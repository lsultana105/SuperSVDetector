use anyhow::Result;
use mpi::traits::*;
use serde::{Serialize, Deserialize};
use bincode;

use crate::{bins, process_bin};
use crate::io_utils;

/// Serializable version of Bin for MPI
#[derive(Serialize, Deserialize)]
struct BinMsg {
    chrom: String,
    start: u64,
    end: u64,
}

impl From<bins::Bin> for BinMsg {
    fn from(b: bins::Bin) -> Self {
        BinMsg { chrom: b.chrom, start: b.start, end: b.end }
    }
}

impl From<BinMsg> for bins::Bin {
    fn from(m: BinMsg) -> Self {
        bins::Bin { chrom: m.chrom, start: m.start, end: m.end }
    }
}

/// Master-slave MPI distribution: master distributes, slaves process bins
pub fn run_mpi(
    ref_fa: &str,
    bam: &str,
    bins_path: &str,
    outdir: &str,
    band: usize,
    k: usize,
    w: usize,
) -> Result<()> {
    let universe = mpi::initialize().expect("MPI init failed");
    let world = universe.world();
    let rank = world.rank();
    let size = world.size();

    if size < 2 {
        panic!("MPI run requires at least 2 ranks (1 master + 1 slave)");
    }

    if rank == 0 {
        // ---- MASTER ----
        let all_bins = io_utils::read_bins(bins_path)?;

        // distribute bins to slaves (ranks 1..size-1)
        for (i, b) in all_bins.into_iter().enumerate() {
            let target = 1 + (i % (size as usize - 1)) as i32;
            let msg: BinMsg = b.into();
            let encoded = bincode::serialize(&msg)?;
            world.process_at_rank(target).send(&encoded[..]);
        }

        // send "done" signal to all slaves
        for target in 1..size {
            world.process_at_rank(target).send(&[] as &[u8]);
        }

        println!("Master finished distributing bins.");
        world.barrier(); // wait for all slaves to finish
    } else {
        // ---- SLAVES ----
        let mut bins_to_process: Vec<bins::Bin> = Vec::new();
        loop {
            let (encoded, _status) = world.process_at_rank(0).receive_vec::<u8>();
            if encoded.is_empty() {
                break; // end signal
            }
            let msg: BinMsg = bincode::deserialize(&encoded)?;
            bins_to_process.push(msg.into());
        }

        println!("Rank {} received {} bins to process", rank, bins_to_process.len());

        // process all assigned bins
        for b in bins_to_process {
            println!("Rank {} processing bin {}:{}-{}", rank, b.chrom, b.start, b.end);
            if let Err(e) = process_bin(ref_fa, bam, outdir, b, band, k, w) {
                eprintln!("Error processing bin on rank {}: {}", rank, e);
            }
        }

        println!("Rank {} finished all assigned bins", rank);
        world.barrier();
    }

    Ok(())
}
