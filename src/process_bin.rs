use anyhow::{Context, Result};

use crate::bins::Bin;
use crate::confirm::{confirm_breakpoint, SvCall};
use crate::hotspot::{hotspot_near_boundary, Hotspot};
use crate::io_utils::{
    fetch_reads_overlapping_with_reader,
    fetch_reference_window_with_reader,
    get_contig_len_from_fai,
    WorkerIo,
};
use crate::kmer::{discover_hotspots, MinimizerIndex};
use crate::suffix::SuffixWorkspace;

#[derive(Clone, Debug)]
pub struct CallConfig {
    pub band: usize,
    pub k: usize,
    pub w: usize,
    pub insert_mean: f64,
    pub insert_sd: f64,
    pub edge_margin: usize,
    pub out_dir: String,
    pub ref_fa: String,
}

pub fn process_bin_with_readers(
    io: &mut WorkerIo,
    bin: &Bin,
    cfg: &CallConfig,
) -> Result<Vec<SvCall>> {
    let chrom = &bin.chrom;

    let debug_enabled = std::env::var("SSV_DEBUG").is_ok();
    let debug_region = std::env::var("SSV_DEBUG_REGION").ok();
    let debug_this_bin = if !debug_enabled {
        false
    } else if let Some(spec) = debug_region.as_deref() {
        let mut ok = false;
        if let Some((chr, rest)) = spec.split_once(':') {
            if chr == chrom {
                if let Some((s, e)) = rest.split_once('-') {
                    if let (Ok(rs), Ok(re)) = (s.parse::<u64>(), e.parse::<u64>()) {
                        ok = bin.start < re && bin.end > rs;
                    }
                }
            }
        }
        ok
    } else {
        true
    };

    // 1) initial fetch
    let ref_window0 = fetch_reference_window_with_reader(
        &mut io.fasta,
        chrom,
        bin.start,
        bin.end,
    )
        .with_context(|| format!("failed reference fetch for {}:{}-{}", chrom, bin.start, bin.end))?;

    let reads0 = fetch_reads_overlapping_with_reader(io, chrom, bin.start, bin.end)
        .with_context(|| format!("failed BAM fetch for {}:{}-{}", chrom, bin.start, bin.end))?;

    if debug_this_bin {
        log::warn!(
            "[BIN] {}:{}-{} reads={} insert_mean={} insert_sd={}",
            chrom, bin.start, bin.end, reads0.len(), cfg.insert_mean, cfg.insert_sd
        );
    }

    let mut sfx = SuffixWorkspace::build(&ref_window0);

    // 2) hotspot discovery
    let mindex0 = MinimizerIndex::build(&ref_window0, cfg.k, cfg.w);
    let mut hotspots: Vec<Hotspot> = discover_hotspots(
        &reads0,
        &mindex0,
        cfg.k,
        bin.start,
        cfg.insert_mean,
        cfg.insert_sd,
    );

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

    // 3) optional boundary expansion
    let needs_expansion = hotspot_near_boundary(&hotspots, ref_window0.len(), cfg.edge_margin);

    let (final_ref_window, final_bin_start, final_reads) = if needs_expansion {
        let contig_len = get_contig_len_from_fai(&cfg.ref_fa, chrom)?;
        let new_start = bin.start.saturating_sub(cfg.edge_margin as u64);
        let new_end = (bin.end + cfg.edge_margin as u64).min(contig_len);

        if debug_this_bin {
            log::warn!(
                "[BIN] expanding: old=[{}-{}] new=[{}-{}]",
                bin.start, bin.end, new_start, new_end
            );
        }

        let extended_ref = fetch_reference_window_with_reader(
            &mut io.fasta,
            chrom,
            new_start,
            new_end,
        )?;

        if new_start < bin.start {
            let shift = (bin.start - new_start) as usize;
            for hs in &mut hotspots {
                hs.local_pos += shift;
            }
        }

        let reads1 = fetch_reads_overlapping_with_reader(io, chrom, new_start, new_end)?;
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

    // 4) confirm hotspots
    let mut calls = Vec::new();

    for hs in &hotspots {
        if let Some(call) = confirm_breakpoint(
            final_bin_start,
            chrom,
            &final_ref_window,
            &sfx,
            hs,
            cfg.band,
            cfg.insert_mean,
            cfg.insert_sd,
            &final_reads,
            cfg.k,
            &io.tid_names,
        )? {
            calls.push(call);
        }
    }

    if debug_this_bin {
        log::warn!("[BIN] calls_written={}", calls.len());
    }

    Ok(calls)
}