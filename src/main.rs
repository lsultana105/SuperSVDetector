use anyhow::Result;
use clap::{Parser, Subcommand};
use log::*;
use rayon::prelude::*;

mod bins;
mod io_utils;
mod kmer;
mod suffix;
mod confirm;
mod mpi_wrap;
mod vcfout;
mod hotspot;
mod cluster;
mod merge_vcf;

use crate::bins::Bin;
use crate::kmer::{discover_hotspots, MinimizerIndex};
use crate::hotspot::{Hotspot, hotspot_near_boundary};
use crate::suffix::SuffixWorkspace;
use crate::confirm::confirm_breakpoint;

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

        // FIXED: Increased k to 31 for better specificity
        #[arg(long, default_value_t = 31)]
        k: usize,

        // FIXED: Increased w to 32 to satisfy the w >= k requirement
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

            if mpi {
                mpi_wrap::run_mpi(&ref_fa, &bam, &bins_path, &outdir, band, k, w)?;
            } else {
                let bins = io_utils::read_bins(&bins_path)?;
                bins.into_par_iter().for_each(|bin| {
                    if let Err(e) =
                        process_bin(&ref_fa, &bam, &outdir, bin, band, k, w)
                    {
                        error!("bin failed: {e}");
                    }
                });
            }
        }

        Command::Merge { outdir, vcf, ref_fa } => {
            vcfout::merge_json_to_vcf(&outdir, &vcf, &ref_fa)?;
            info!("VCF written to {vcf}");
        }
    }

    Ok(())
}

pub fn process_bin(
    ref_fa: &str,
    bam_path: &str,
    outdir: &str,
    mut bin: Bin,
    band: usize,
    k: usize,
    w: usize,
) -> Result<()> {
    use crate::io_utils::*;

    const EDGE_MARGIN: usize = 2000;

    let ref_window = fetch_reference_window(ref_fa, &bin.chrom, bin.start, bin.end)?;
    let sfx = SuffixWorkspace::build(&ref_window);

    let mut reader = open_bam(bam_path)?;
    let reads = fetch_reads_overlapping(&mut reader, &bin.chrom, bin.start, bin.end)?;

    let mindex = MinimizerIndex::build(&ref_window, k, w);
    let mut hotspots = discover_hotspots(&reads, &mindex, k);

    let needs_expansion =
        hotspot_near_boundary(&hotspots[..], ref_window.len(), EDGE_MARGIN);

    let (final_ref_window, final_bin_start) = if needs_expansion {
        let new_start = bin.start.saturating_sub(EDGE_MARGIN as u64);
        let new_end = bin.end + EDGE_MARGIN as u64;

        let extended = fetch_reference_window(ref_fa, &bin.chrom, new_start, new_end)?;

        for hs in &mut hotspots {
            hs.local_pos += (bin.start - new_start) as usize;
        }

        (extended, new_start)
    } else {
        (ref_window.clone(), bin.start)
    };

    let calls: Vec<_> = hotspots
        .into_par_iter()
        .filter_map(|hs| {
            confirm_breakpoint(
                final_bin_start,
                &bin.chrom,
                &final_ref_window,
                &sfx,
                &hs,
                band,
            )
                .ok()
                .flatten()
        })
        .collect();

    vcfout::write_bin_json(outdir, &bin, &calls)?;
    Ok(())
}

