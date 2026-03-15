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

fn collapse_nearby_hotspots(mut hotspots: Vec<Hotspot>) -> Vec<Hotspot> {
    if hotspots.is_empty() {
        return hotspots;
    }

    hotspots.sort_by(|a, b| {
        let ord = a.reason.cmp(&b.reason);
        if ord == std::cmp::Ordering::Equal {
            a.local_pos.cmp(&b.local_pos)
        } else {
            ord
        }
    });

    let mut out = Vec::new();
    let mut cur = hotspots[0].clone();

    const HOTSPOT_MERGE_WIN: usize = 300;

    for hs in hotspots.into_iter().skip(1) {
        if hs.reason == cur.reason && hs.local_pos.abs_diff(cur.local_pos) <= HOTSPOT_MERGE_WIN {
            cur.local_pos = (cur.local_pos + hs.local_pos) / 2;
        } else {
            out.push(cur);
            cur = hs;
        }
    }

    out.push(cur);
    out.sort_by_key(|h| h.local_pos);
    out
}

fn same_call(a: &SvCall, b: &SvCall) -> bool {
    if a.chrom != b.chrom || a.svtype != b.svtype {
        return false;
    }

    match a.svtype.as_str() {
        "DEL" => a.pos.abs_diff(b.pos) <= 100 && a.end.abs_diff(b.end) <= 100,
        "INS" => a.pos.abs_diff(b.pos) <= 50 && (a.len - b.len).abs() <= 50,
        "BND" => {
            a.pos.abs_diff(b.pos) <= 100
                && a.mate_chrom == b.mate_chrom
                && match (a.mate_pos, b.mate_pos) {
                (Some(x), Some(y)) => x.abs_diff(y) <= 100,
                (None, None) => true,
                _ => false,
            }
        }
        _ => a.pos.abs_diff(b.pos) <= 100 && a.end.abs_diff(b.end) <= 100,
    }
}

fn dedup_calls(mut calls: Vec<SvCall>) -> Vec<SvCall> {
    if calls.is_empty() {
        return calls;
    }

    calls.sort_by(|a, b| {
        let ord = a.chrom.cmp(&b.chrom);
        if ord != std::cmp::Ordering::Equal {
            return ord;
        }
        let ord = a.svtype.cmp(&b.svtype);
        if ord != std::cmp::Ordering::Equal {
            return ord;
        }
        let ord = a.pos.cmp(&b.pos);
        if ord != std::cmp::Ordering::Equal {
            return ord;
        }
        a.end.cmp(&b.end)
    });

    let mut out: Vec<SvCall> = Vec::new();

    for call in calls {
        if let Some(last) = out.last_mut() {
            if same_call(last, &call) {
                if call.score > last.score {
                    *last = call;
                }
                continue;
            }
        }
        out.push(call);
    }

    out
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

    // 1) Initial fetch
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
            chrom,
            bin.start,
            bin.end,
            reads0.len(),
            cfg.insert_mean,
            cfg.insert_sd
        );
    }

    // 2) Discover hotspots once on the original bin
    let mindex0 = MinimizerIndex::build(&ref_window0, cfg.k, cfg.w);
    let hotspots0 = discover_hotspots(
        &reads0,
        &mindex0,
        cfg.k,
        bin.start,
        cfg.insert_mean,
        cfg.insert_sd,
    );
    let mut hotspots = collapse_nearby_hotspots(hotspots0);

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
            "[BIN] hotspots after collapse={} discordant_del={} bnd_pe={} soft_clip_sr={}",
            hotspots.len(),
            d,
            b,
            s
        );
    }

    // 3) Decide whether expansion is needed
    let needs_expansion = hotspot_near_boundary(&hotspots, ref_window0.len(), cfg.edge_margin);

    // 4) Prepare final reference / reads for confirmation
    let (final_ref_window, final_bin_start, final_reads) = if needs_expansion {
        let contig_len = get_contig_len_from_fai(&cfg.ref_fa, chrom)?;
        let new_start = bin.start.saturating_sub(cfg.edge_margin as u64);
        let new_end = (bin.end + cfg.edge_margin as u64).min(contig_len);

        if debug_this_bin {
            log::warn!(
                "[BIN] expanding: old=[{}-{}] new=[{}-{}]",
                bin.start,
                bin.end,
                new_start,
                new_end
            );
        }

        let extended_ref = fetch_reference_window_with_reader(
            &mut io.fasta,
            chrom,
            new_start,
            new_end,
        )
            .with_context(|| format!("failed expanded reference fetch for {}:{}-{}", chrom, new_start, new_end))?;

        if new_start < bin.start {
            let shift = (bin.start - new_start) as usize;
            for hs in &mut hotspots {
                hs.local_pos += shift;
            }
        }

        let reads1 = fetch_reads_overlapping_with_reader(io, chrom, new_start, new_end)
            .with_context(|| format!("failed expanded BAM fetch for {}:{}-{}", chrom, new_start, new_end))?;

        if debug_this_bin {
            log::warn!(
                "[BIN] expanded ref_len={} reads={} hotspots_preserved={}",
                extended_ref.len(),
                reads1.len(),
                hotspots.len()
            );
        }

        (extended_ref, new_start, reads1)
    } else {
        (ref_window0, bin.start, reads0)
    };

    let sfx = SuffixWorkspace::build(&final_ref_window);

    // 5) Confirm hotspots
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

    let calls = dedup_calls(calls);

    if debug_this_bin {
        log::warn!("[BIN] calls_after_dedup={}", calls.len());
    }

    Ok(calls)
}