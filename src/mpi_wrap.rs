use anyhow::Result;
use mpi::traits::*;
use serde::{Deserialize, Serialize};

use crate::{bins, io_utils};

#[derive(Serialize, Deserialize)]
struct BinMsg {
    chrom: String,
    start: u64,
    end: u64,
}

impl From<bins::Bin> for BinMsg {
    fn from(b: bins::Bin) -> Self {
        BinMsg {
            chrom: b.chrom,
            start: b.start,
            end: b.end,
        }
    }
}

pub fn run_mpi(
    ref_fa: &str,
    bam: &str,
    bins_path: &str,
    outdir: &str,
    band: usize,
    k: usize,
    w: usize,
    insert_mean: f64,
    insert_sd: f64,
    tid_names: &[String],
) -> Result<()> {
    let universe = mpi::initialize().expect("MPI init failed");
    let world = universe.world();
    let (rank, size) = (world.rank(), world.size());

    // ✅ FIX: if size==1, run locally on rank 0 (no workers exist)
    if size == 1 {
        if rank == 0 {
            eprintln!("[RANK 0] size==1: running bins locally (no MPI workers).");
            let bins = io_utils::read_bins(bins_path)?;
            for bin in bins {
                if let Err(e) = crate::process_bin(
                    ref_fa,
                    bam,
                    outdir,
                    bin,
                    band,
                    k,
                    w,
                    insert_mean,
                    insert_sd,
                    tid_names,
                ) {
                    eprintln!("[RANK 0] bin failed: {}", e);
                }
            }
        }
        return Ok(());
    }

    // Master
    if rank == 0 {
        eprintln!("[MASTER] Initializing... Total Ranks: {}", size);
        eprintln!("[MASTER] Reading bins from: {}", bins_path);

        let mut all_bins = io_utils::read_bins(bins_path)?;
        eprintln!("[MASTER] Loaded {} bins. Starting distribution...", all_bins.len());

        let mut finished_workers = 0;

        while finished_workers < (size - 1) {
            // workers send "ready" pings
            let (_, status) = world.any_process().receive_vec::<u8>();
            let worker = status.source_rank();

            if let Some(bin) = all_bins.pop() {
                let msg = BinMsg::from(bin);
                let encoded = bincode::serialize(&msg)?;
                world.process_at_rank(worker).send(&encoded[..]);
            } else {
                eprintln!("[MASTER] No more bins. Sending shutdown to Rank {}", worker);
                world.process_at_rank(worker).send(&[] as &[u8]);
                finished_workers += 1;
            }
        }

        eprintln!("[MASTER] All tasks complete. Finalizing...");
        return Ok(());
    }

    // Workers
    eprintln!("[RANK {}] Worker online.", rank);
    loop {
        // ping master that we are ready
        world.process_at_rank(0).send(&[1u8] as &[u8]);

        let (encoded, _) = world.process_at_rank(0).receive_vec::<u8>();
        if encoded.is_empty() {
            eprintln!("[RANK {}] Shutting down.", rank);
            break;
        }

        let msg: BinMsg = bincode::deserialize(&encoded)?;
        let bin = bins::Bin {
            chrom: msg.chrom,
            start: msg.start,
            end: msg.end,
        };

        eprintln!(
            "[RANK {}] Working on {}:{}-{}",
            rank, bin.chrom, bin.start, bin.end
        );

        if let Err(e) = crate::process_bin(
            ref_fa,
            bam,
            outdir,
            bin,
            band,
            k,
            w,
            insert_mean,
            insert_sd,
            tid_names,
        ) {
            eprintln!("[RANK {}] ERROR in process_bin: {}", rank, e);
        }
    }

    Ok(())
}
