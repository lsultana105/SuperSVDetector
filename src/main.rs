use anyhow::Result;
use clap::{Parser, Subcommand};
use log::*;
use rayon::prelude::*;
use std::time::Instant;
use crate::io_utils::WorkerIo;
use crate::process_bin::{process_bin, CallConfig};

mod bins;
mod io_utils;
mod kmer;
mod suffix;
mod confirm;
mod mpi_wrap;
mod vcfout;
mod hotspot;
mod cluster;
mod process_bin;

// The pipeline is split into three stages:
// 1. index reference bins
// 2. call SV candidates per bin
// 3. merge per-bin calls into a final VCF
#[derive(Parser, Debug)]
#[command(
    name = "SuperSVDetector",
    about = "Structural Variant Detection using minimizers and suffix arrays"
)]
struct Cli {
    #[command(subcommand)]
    cmd: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    Index {
        #[arg(long)]
        ref_fa: String,

        #[arg(long)]
        bin: usize,

        #[arg(long, default_value_t = 0)]
        pad: usize,

        #[arg(long)]
        out: String,
    },

    Call {
        #[arg(long)]
        ref_fa: String,

        #[arg(long)]
        bam: String,

        #[arg(long)]
        bins: String,

        #[arg(long)]
        outdir: String,

        #[arg(long, default_value_t = 100)]
        band: usize,

        #[arg(long, default_value_t = 31)]
        k: usize,

        #[arg(long, default_value_t = 32)]
        w: usize,

        #[arg(long)]
        mpi: bool,

        #[arg(long)]
        threads: Option<usize>,
    },

    Merge {
        #[arg(long)]
        outdir: String,

        #[arg(long)]
        vcf: String,

        #[arg(long)]
        ref_fa: String,
    },
}

fn main() -> Result<()> {
    let start = Instant::now();

    env_logger::init();
    let cli = Cli::parse();

    match cli.cmd {
        Command::Index { ref_fa, bin, pad, out } => {
            bins::build_bins(&ref_fa, bin, pad, &out)?;
            info!("Bins written to {out}");
        }

        Command::Call {
            ref_fa,
            bam,
            bins: bins_path,
            outdir,
            band,
            k,
            w,
            mpi,
            threads,
        } => {
            if let Some(t) = threads {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(t)
                    .build_global()
                    .ok();
            }

            let (insert_mean, insert_sd) =
                io_utils::estimate_insert_stats(&bam).unwrap_or((0.0, 0.0));
            info!("Insert stats mean={:.2} sd={:.2}", insert_mean, insert_sd);

            let bins = io_utils::read_bins(&bins_path)?;

            let cfg = CallConfig {
                band,
                k,
                w,
                insert_mean,
                insert_sd,
                edge_margin: 2000,
                out_dir: outdir.clone(),
                ref_fa: ref_fa.clone(),
            };

            if mpi {
                mpi_wrap::run_mpi(&bam, &ref_fa, bins, cfg)?;
            } else {
                bins.into_par_iter().for_each_init(
                    || WorkerIo::new(&bam, &ref_fa).expect("failed to initialize WorkerIo"),
                    |io, bin| {
                        match process_bin(io, &bin, &cfg) {
                            Ok(calls) => {
                                if let Err(e) = vcfout::write_bin_json(&outdir, &bin, &calls) {
                                    error!(
                                        "failed to write JSON for {}:{}-{}: {e}",
                                        bin.chrom, bin.start, bin.end
                                    );
                                }
                            }
                            Err(e) => {
                                error!(
                                    "bin failed for {}:{}-{}: {e}",
                                    bin.chrom, bin.start, bin.end
                                );
                            }
                        }
                    },
                );
            }
        }

        Command::Merge { outdir, vcf, ref_fa } => {
            vcfout::merge_json_to_vcf(&outdir, &vcf, &ref_fa)?;
            info!("VCF written to {vcf}");
        }
    }

    let duration = start.elapsed();
    println!("--- Analysis Complete ---");
    println!("Total Time Taken: {:?}", duration);

    Ok(())
}