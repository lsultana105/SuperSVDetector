use anyhow::{anyhow, Result};
use mpi::traits::*;
use serde::{Deserialize, Serialize};

use crate::bins::Bin;
use crate::io_utils::WorkerIo;
use crate::process_bin::{process_bin_with_readers, CallConfig};
use crate::vcfout;

const TAG_WORK: i32 = 1;
const TAG_STOP: i32 = 2;
const TAG_DONE: i32 = 3;
const TAG_ERR: i32 = 4;

#[derive(Serialize, Deserialize, Debug)]
struct WorkerDone {
    rank: i32,
    chrom: String,
    start: u64,
    end: u64,
}

#[derive(Serialize, Deserialize, Debug)]
struct WorkerErr {
    rank: i32,
    chrom: String,
    start: u64,
    end: u64,
    message: String,
}

pub fn run_mpi(
    bam_path: &str,
    fasta_path: &str,
    bins: Vec<Bin>,
    cfg: CallConfig,
) -> Result<()> {
    let universe = mpi::initialize().ok_or_else(|| anyhow!("failed to initialize MPI"))?;
    let world = universe.world();
    let rank = world.rank();
    let size = world.size();

    if size < 2 {
        return Err(anyhow!(
            "MPI run requires at least 2 ranks (1 master + 1 worker)"
        ));
    }

    if rank == 0 {
        run_master(&world, bins)
    } else {
        run_worker(&world, bam_path, fasta_path, cfg)
    }
}

fn run_master<C: mpi::topology::Communicator>(world: &C, bins: Vec<Bin>) -> Result<()> {
    let worker_count = world.size() - 1;
    let total = bins.len();
    let mut next_idx = 0usize;

    for worker_rank in 1..=worker_count {
        if next_idx < total {
            send_bin(world, worker_rank, &bins[next_idx])?;
            next_idx += 1;
        } else {
            send_stop(world, worker_rank);
        }
    }

    let mut completed = 0usize;

    while completed < total {
        let (msg, status) = world.any_process().receive_vec::<u8>();
        let src = status.source_rank();
        let tag = status.tag();

        match tag {
            TAG_DONE => {
                let done: WorkerDone = serde_json::from_slice(&msg)?;
                eprintln!(
                    "[master] done from rank {} -> {}:{}-{}",
                    done.rank, done.chrom, done.start, done.end
                );

                completed += 1;

                if next_idx < total {
                    send_bin(world, src, &bins[next_idx])?;
                    next_idx += 1;
                } else {
                    send_stop(world, src);
                }
            }
            TAG_ERR => {
                let err: WorkerErr = serde_json::from_slice(&msg)?;
                eprintln!(
                    "[master] ERROR from rank {} on {}:{}-{} -> {}",
                    err.rank, err.chrom, err.start, err.end, err.message
                );

                completed += 1;

                if next_idx < total {
                    send_bin(world, src, &bins[next_idx])?;
                    next_idx += 1;
                } else {
                    send_stop(world, src);
                }
            }
            _ => {
                return Err(anyhow!("master received unexpected MPI tag {}", tag));
            }
        }
    }

    Ok(())
}

fn run_worker<C: mpi::topology::Communicator>(
    world: &C,
    bam_path: &str,
    fasta_path: &str,
    cfg: CallConfig,
) -> Result<()> {
    let rank = world.rank();

    // Open BAM + FASTA once per worker process
    let mut io = WorkerIo::new(bam_path, fasta_path)?;

    loop {
        let (msg, status) = world.process_at_rank(0).receive_vec::<u8>();
        let tag = status.tag();

        match tag {
            TAG_WORK => {
                let bin: Bin = serde_json::from_slice(&msg)?;

                match process_bin_with_readers(&mut io, &bin, &cfg) {
                    Ok(calls) => {
                        vcfout::write_bin_json(&cfg.out_dir, &bin, &calls)?;

                        let done = WorkerDone {
                            rank,
                            chrom: bin.chrom.clone(),
                            start: bin.start,
                            end: bin.end,
                        };
                        let payload = serde_json::to_vec(&done)?;
                        world.process_at_rank(0).send_with_tag(&payload[..], TAG_DONE);
                    }
                    Err(e) => {
                        let err = WorkerErr {
                            rank,
                            chrom: bin.chrom.clone(),
                            start: bin.start,
                            end: bin.end,
                            message: e.to_string(),
                        };
                        let payload = serde_json::to_vec(&err)?;
                        world.process_at_rank(0).send_with_tag(&payload[..], TAG_ERR);
                    }
                }
            }
            TAG_STOP => break,
            _ => {
                return Err(anyhow!(
                    "worker rank {} received unexpected MPI tag {}",
                    rank,
                    tag
                ));
            }
        }
    }

    Ok(())
}

fn send_bin<C: mpi::topology::Communicator>(
    world: &C,
    worker_rank: i32,
    bin: &Bin,
) -> Result<()> {
    let payload = serde_json::to_vec(bin)?;
    world
        .process_at_rank(worker_rank)
        .send_with_tag(&payload[..], TAG_WORK);
    Ok(())
}

fn send_stop<C: mpi::topology::Communicator>(world: &C, worker_rank: i32) {
    let empty: [u8; 0] = [];
    world
        .process_at_rank(worker_rank)
        .send_with_tag(&empty[..], TAG_STOP);
}