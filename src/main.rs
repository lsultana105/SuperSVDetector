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
use std::time::Instant;

use crate::io_utils::fetch_reference_window;
use crate::io_utils::open_bam;
use crate::io_utils::fetch_reads_overlapping;
use rust_htslib::bam::Read;

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

            let (insert_mean, insert_sd) = io_utils::estimate_insert_stats(&bam).unwrap_or((0.0, 0.0));
            info!("Insert stats mean={:.2} sd={:.2}", insert_mean, insert_sd);

            let mut reader = rust_htslib::bam::IndexedReader::from_path(&bam)?;
            let header = reader.header().to_owned();

            let mut tid_names: Vec<String> = Vec::new();
            for tid in 0..header.target_count() {
                let name = String::from_utf8_lossy(header.tid2name(tid)).to_string();
                tid_names.push(name);
            }

            if mpi {
                mpi_wrap::run_mpi(
                    &ref_fa,
                    &bam,
                    &bins_path,
                    &outdir,
                    band,
                    k,
                    w,
                    insert_mean,
                    insert_sd,
                    &tid_names,
                )?;
            } else {
                let mut reader = rust_htslib::bam::IndexedReader::from_path(&bam)?;
                let header = reader.header().to_owned();

                let mut tid_names: Vec<String> = Vec::new();
                for tid in 0..header.target_count() {
                    let name = String::from_utf8_lossy(header.tid2name(tid)).to_string();
                    tid_names.push(name);
                }

                let bins = io_utils::read_bins(&bins_path)?;
                bins.into_par_iter().for_each(|bin| {
                    if let Err(e) =
                        process_bin(&ref_fa, &bam, &outdir, bin, band, k, w, insert_mean, insert_sd, &tid_names,)
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
    let duration = start.elapsed();
    println!("--- Analysis Complete ---");
    println!("Total Time Taken: {:?}", duration);

    Ok(())
}

pub fn process_bin(
    ref_fa: &str,
    bam_path: &str,
    outdir: &str,
    bin: Bin,
    band: usize,
    k: usize,
    w: usize,
    insert_mean: f64,
    insert_sd: f64,
    tid_names: &[String],
) -> Result<()> {
    use crate::io_utils::*;
    use rayon::prelude::*;

    const EDGE_MARGIN: usize = 2000;

    // -------- Optional generic debug (NO hardcoding) --------
    let debug_enabled = std::env::var("SSV_DEBUG").is_ok();

    // If set, only debug bins overlapping this region: "chr:start-end"
    let debug_region = std::env::var("SSV_DEBUG_REGION").ok();
    let debug_this_bin = if !debug_enabled {
        false
    } else if let Some(spec) = debug_region.as_deref() {
        // parse "chr1:119000000-120000000"
        let mut ok = false;
        if let Some((chr, rest)) = spec.split_once(':') {
            if chr == bin.chrom {
                if let Some((s, e)) = rest.split_once('-') {
                    if let (Ok(rs), Ok(re)) = (s.parse::<u64>(), e.parse::<u64>()) {
                        // overlap check
                        ok = bin.start < re && bin.end > rs;
                    }
                }
            }
        }
        ok
    } else {
        true
    };
    // --------------------------------------------------------

    // 1) Fetch reference for original bin
    let ref_window0 = fetch_reference_window(ref_fa, &bin.chrom, bin.start, bin.end)?;
    let mut sfx = SuffixWorkspace::build(&ref_window0);

    // 2) Fetch reads for original bin
    let mut reader0 = open_bam(bam_path)?;
    let reads0 = fetch_reads_overlapping(&mut reader0, &bin.chrom, bin.start, bin.end)?;

    if debug_this_bin {
        log::warn!(
            "[BIN] {}:{}-{} reads={} insert_mean={} insert_sd={}",
            bin.chrom, bin.start, bin.end, reads0.len(), insert_mean, insert_sd
        );
    }

    // 3) Discover hotspots (on original window + original reads)
    let mindex0 = MinimizerIndex::build(&ref_window0, k, w);
    let mut hotspots = discover_hotspots(&reads0, &mindex0, k, bin.start, insert_mean, insert_sd);

    if debug_this_bin {
        let mut d = 0usize;
        let mut b = 0usize;
        let mut s = 0usize;
        for h in &hotspots {
            match h.reason.as_str() {
                "discordant_del" => d += 1,
                "bnd_pe" => b += 1,
                "soft_clip_sr" => s += 1,
                _ => {}
            }
        }
        log::warn!(
            "[BIN] hotspots total={} discordant_del={} bnd_pe={} soft_clip_sr={}",
            hotspots.len(),
            d,
            b,
            s
        );
    }

    // 4) Expand if any hotspot near boundary
    let needs_expansion = hotspot_near_boundary(&hotspots[..], ref_window0.len(), EDGE_MARGIN);

    // Finalize (ref_window, bin_start, reads) so they match each other.
    let (final_ref_window, final_bin_start, final_reads) = if needs_expansion {
        let new_start = bin.start.saturating_sub(EDGE_MARGIN as u64);
        let new_end = bin.end + EDGE_MARGIN as u64;

        if debug_this_bin {
            log::warn!(
                "[BIN] expanding: old=[{}-{}] new=[{}-{}]",
                bin.start, bin.end, new_start, new_end
            );
        }

        let extended_ref = fetch_reference_window(ref_fa, &bin.chrom, new_start, new_end)?;

        // shift hotspot local_pos because local coordinate system changed
        if new_start < bin.start {
            let shift = (bin.start - new_start) as usize;
            for hs in &mut hotspots {
                hs.local_pos += shift;
            }
        }

        // ✅ CRITICAL: refetch reads for expanded region too
        let mut reader1 = open_bam(bam_path)?;
        let reads1 = fetch_reads_overlapping(&mut reader1, &bin.chrom, new_start, new_end)?;

        // suffix must match expanded ref
        sfx = SuffixWorkspace::build(&extended_ref);

        if debug_this_bin {
            log::warn!(
                "[BIN] expanded ref_len={} reads={}",
                extended_ref.len(),
                reads1.len()
            );
        }

        (extended_ref, new_start, reads1)
    } else {
        (ref_window0, bin.start, reads0)
    };

    // 5) Confirm hotspots in parallel
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
                insert_mean,
                insert_sd,
                &final_reads,
                k,
                tid_names,
            )
                .ok()
                .flatten()
        })
        .collect();

    if debug_this_bin {
        log::warn!("[BIN] calls_written={}", calls.len());
    }

    vcfout::write_bin_json(outdir, &bin, &calls)?;
    Ok(())
}




