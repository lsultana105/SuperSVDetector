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

/// Simple MPI distribution for demo:
/// - Root reads all bins
/// - Each rank gets bins in round-robin
/// - Each rank processes its bins
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

    // Root reads bins
    let all_bins = if rank == 0 {
        io_utils::read_bins(bins_path)?
    } else {
        Vec::new()
    };

    // Assign bins
    let bins_to_process: Vec<bins::Bin> = if rank == 0 {
        let mut assignments = vec![Vec::new(); size as usize];
        for (i, b) in all_bins.into_iter().enumerate() {
            let target = i % size as usize;
            if target == 0 {
                assignments[0].push(b);
            } else {
                let msg: BinMsg = b.into();
                let encoded = bincode::serialize(&msg).unwrap();
                world.process_at_rank(target as i32).send(&encoded[..]);
            }
        }

        // Send an empty Vec<u8> to signal "done"
        for target in 1..size {
            world.process_at_rank(target as i32).send(&[] as &[u8]);
        }

        assignments[0].clone()
    } else {
        let mut received = Vec::new();
        loop {
            let (encoded, _status) = world.process_at_rank(0).receive_vec::<u8>();

            if encoded.is_empty() {
                break; // end signal
            }

            let msg: BinMsg = bincode::deserialize(&encoded).unwrap();
            received.push(msg.into());
        }
        received
    };

    println!("MPI rank {} started with {} bins", rank, bins_to_process.len());

    // Process bins
    for b in bins_to_process {
        println!("Rank {} processing bin {}:{}-{}", rank, b.chrom, b.start, b.end);
        if let Err(e) = process_bin(ref_fa, bam, outdir, b, band, k, w) {
            eprintln!("Error processing bin: {}", e);
        }
    }

    println!("Rank {} finished all assigned bins", rank);
    world.barrier();
    Ok(())
}
