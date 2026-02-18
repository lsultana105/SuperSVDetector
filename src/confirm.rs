use anyhow::Result;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

use crate::hotspot::Hotspot;
use crate::suffix::SuffixWorkspace;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SvCall {
    pub chrom: String,
    pub pos: u64,
    pub end: u64,
    pub svtype: String,
    pub len: isize,
    pub score: i32,
    pub reason: String,

    // Minimal BND support (optional)
    pub mate_chrom: Option<String>,
    pub mate_pos: Option<u64>,
}

// ------------------------
// small utilities
// ------------------------

fn median_u64(v: &mut Vec<u64>) -> u64 {
    v.sort_unstable();
    v[v.len() / 2]
}

fn median_isize(v: &mut Vec<isize>) -> isize {
    v.sort_unstable();
    v[v.len() / 2]
}

fn ref_consumed_len(cig: &bam::record::CigarStringView) -> u64 {
    use bam::record::Cigar::*;
    cig.iter()
        .map(|c| match c {
            Match(n) | Diff(n) | Equal(n) | Del(n) | RefSkip(n) => *n as u64,
            Ins(_) | SoftClip(_) | HardClip(_) | Pad(_) => 0,
        })
        .sum()
}

/// Spread after trimming extremes (trim_frac on each side).
/// If n is small, returns 0 (i.e. don't over-penalize small samples).
fn trimmed_spread(sorted: &[u64], trim_frac: f64) -> u64 {
    let n = sorted.len();
    if n < 5 {
        return 0;
    }
    let lo = (n as f64 * trim_frac).floor() as usize;
    let mut hi = (n as f64 * (1.0 - trim_frac)).ceil() as usize;
    if hi >= n {
        hi = n - 1;
    }
    if hi <= lo {
        return 0;
    }
    sorted[hi] - sorted[lo]
}

// ------------------------
// PE deletion confirmer (the “good” one)
// ------------------------

fn confirm_pe_deletion_from_support(
    support: &[&bam::Record],
    insert_mean: f64,
    insert_sd: f64,
) -> Option<(u64, u64, usize)> {
    use std::collections::HashSet;

    // Tunables
    const MIN_MAPQ: u8 = 20;
    const MIN_DEL_BP: u64 = 50;

    // keep this at 4 for TBX sensitivity (you already moved toward this)
    const MIN_PE_SUPPORT: usize = 4;

    // keep these aligned with discover_hotspots (you already did Z=3 there in TBX debug)
    const Z_CUTOFF: f64 = 3.0;

    const MAX_DEL_SPAN: u64 = 50_000;

    // robust gates
    const MAX_BP_SPREAD: u64 = 1500;
    const TRIM_FRAC: f64 = 0.10;

    // NEW: breakpoint clustering window (handles bimodal mate starts)
    const CLUSTER_WIN_BP: u64 = 800;

    if insert_sd <= 0.0 {
        return None;
    }

    // Two-pointer densest window on sorted positions (returns [lo, hi) indices)
    fn densest_window(sorted: &[u64], width: u64) -> (usize, usize) {
        if sorted.is_empty() {
            return (0, 0);
        }
        let mut best_lo = 0usize;
        let mut best_hi = 1usize;

        let mut lo = 0usize;
        for hi in 0..sorted.len() {
            while sorted[hi].saturating_sub(sorted[lo]) > width {
                lo += 1;
            }
            if (hi + 1) - lo > best_hi - best_lo {
                best_lo = lo;
                best_hi = hi + 1;
            }
        }
        (best_lo, best_hi)
    }

    let hard_tlen_min: f64 = (insert_mean + 3.0 * insert_sd).max(800.0);

    let mut left_bps: Vec<u64> = Vec::new();
    let mut right_bps: Vec<u64> = Vec::new();
    let mut seen_qname: HashSet<Vec<u8>> = HashSet::new();

    for r in support {
        if r.is_secondary() || r.is_duplicate() {
            continue;
        }
        if r.mapq() < MIN_MAPQ {
            continue;
        }
        if !(r.is_paired() && !r.is_unmapped() && !r.is_mate_unmapped()) {
            continue;
        }
        if r.tid() != r.mtid() {
            continue;
        }

        // discordant test consistent with discover_hotspots: z-pass OR hard cutoff
        let tlen = r.insert_size().abs() as f64;
        let mut discordant = false;

        let z = (tlen - insert_mean) / insert_sd;
        if z >= Z_CUTOFF {
            discordant = true;
        }
        if tlen >= hard_tlen_min {
            discordant = true;
        }
        if !discordant {
            continue;
        }

        // orientation filter (same as you used before)
        let rev = r.is_reverse();
        let mrev = r.is_mate_reverse();
        if rev == mrev {
            continue;
        }

        // unique fragment by qname
        let qn = r.qname().to_vec();
        if !seen_qname.insert(qn) {
            continue;
        }

        let rpos0 = r.pos().max(0) as u64;
        let mpos0 = r.mpos().max(0) as u64;

        // only use the left read in the pair
        if rpos0 > mpos0 {
            continue;
        }

        let left_bp0 = rpos0 + ref_consumed_len(&r.cigar()); // end of left read
        let right_bp0 = mpos0; // start of right read

        if right_bp0 <= left_bp0 {
            continue;
        }
        if (right_bp0 - left_bp0) > MAX_DEL_SPAN {
            continue;
        }

        left_bps.push(left_bp0);
        right_bps.push(right_bp0);
    }

    if left_bps.len() < MIN_PE_SUPPORT || right_bps.len() < MIN_PE_SUPPORT {
        return None;
    }

    left_bps.sort_unstable();
    right_bps.sort_unstable();

    // NEW: pick densest cluster to kill bimodality (TBX problem)
    let (l_lo, l_hi) = densest_window(&left_bps, CLUSTER_WIN_BP);
    let (r_lo, r_hi) = densest_window(&right_bps, CLUSTER_WIN_BP);

    let mut left_c = left_bps[l_lo..l_hi].to_vec();
    let mut right_c = right_bps[r_lo..r_hi].to_vec();

    if left_c.len() < MIN_PE_SUPPORT || right_c.len() < MIN_PE_SUPPORT {
        return None;
    }

    left_c.sort_unstable();
    right_c.sort_unstable();

    // robust spread gate (trimmed)
    let left_spread = trimmed_spread(&left_c, TRIM_FRAC);
    let right_spread = trimmed_spread(&right_c, TRIM_FRAC);
    if left_spread > MAX_BP_SPREAD || right_spread > MAX_BP_SPREAD {
        return None;
    }

    let start0 = median_u64(&mut left_c);
    let end0 = median_u64(&mut right_c);

    if end0 <= start0 || (end0 - start0) < MIN_DEL_BP {
        return None;
    }

    // score = cluster support used (min of both sides)
    let score = left_c.len().min(right_c.len());
    Some((start0, end0, score))
}



// ------------------------
// main entry point
// ------------------------

pub fn confirm_breakpoint(
    bin_start: u64,
    chrom: &str,
    _ref_window: &[u8],
    _sfx: &SuffixWorkspace,
    hs: &Hotspot,
    band: usize,
    insert_mean: f64,
    insert_sd: f64,
    reads: &[bam::Record],
    k: usize,
    tid_names: &[String],
) -> Result<Option<SvCall>> {
    // Window around hotspot (global coords)
    let center = hs.local_pos;
    let local_start = center.saturating_sub(band);
    let local_end = center + band;

    let g_start0 = bin_start + local_start as u64;
    let g_end0 = bin_start + local_end as u64;

    // Gather overlapping reads (or mate overlaps) into support
    let mut support: Vec<&bam::Record> = Vec::new();
    for r in reads.iter() {
        let pos0 = r.pos().max(0) as u64;
        let end0 = pos0 + ref_consumed_len(&r.cigar());

        let read_overlaps = end0 > g_start0 && pos0 < g_end0;

        let mate_overlaps = r.is_paired()
            && !r.is_mate_unmapped()
            && (r.mpos() as u64) >= g_start0
            && (r.mpos() as u64) < g_end0;

        if read_overlaps || mate_overlaps {
            support.push(r);
        }
    }

    if support.is_empty() {
        return Ok(None);
    }

    // -------- DEL from PE --------
    if hs.reason == "discordant_del" {
        if let Some((start0, end0, score)) =
            confirm_pe_deletion_from_support(&support, insert_mean, insert_sd)
        {
            let svlen = -((end0 - start0) as isize);

            return Ok(Some(SvCall {
                chrom: chrom.into(),
                pos: start0 + 1, // VCF 1-based
                end: end0 + 1,
                svtype: "DEL".into(),
                len: svlen,
                score: score as i32,
                reason: "PE".into(),
                mate_chrom: None,
                mate_pos: None,
            }));
        }
        return Ok(None);
    }

    // -------- BND minimal --------
    if hs.reason == "bnd_pe" {
        const MIN_MAPQ: u8 = 20;
        const MIN_BND_SUPPORT: usize = 3;
        const MAX_POS_SPREAD: u64 = 500;

        let mut a_pos: Vec<u64> = Vec::new();
        let mut b_pos: Vec<u64> = Vec::new();
        let mut mate_chrom_name: Option<String> = None;
        let mut seen_qname: HashSet<Vec<u8>> = HashSet::new();

        for r in &support {
            if r.is_secondary() || r.is_duplicate() {
                continue;
            }
            if r.mapq() < MIN_MAPQ {
                continue;
            }
            if !(r.is_paired() && !r.is_unmapped() && !r.is_mate_unmapped()) {
                continue;
            }
            if r.tid() == r.mtid() {
                continue;
            }

            let qn = r.qname().to_vec();
            if !seen_qname.insert(qn) {
                continue;
            }

            let rpos0 = r.pos().max(0) as u64;
            let mpos0 = r.mpos().max(0) as u64;
            a_pos.push(rpos0);
            b_pos.push(mpos0);

            let mtid = r.mtid();
            if mtid >= 0 {
                let mtid = mtid as usize;
                if mtid < tid_names.len() {
                    mate_chrom_name = Some(tid_names[mtid].clone());
                }
            }
        }

        if a_pos.len() < MIN_BND_SUPPORT || mate_chrom_name.is_none() {
            return Ok(None);
        }

        a_pos.sort_unstable();
        b_pos.sort_unstable();

        let a_spread = a_pos[a_pos.len() - 1] - a_pos[0];
        let b_spread = b_pos[b_pos.len() - 1] - b_pos[0];
        if a_spread > MAX_POS_SPREAD || b_spread > MAX_POS_SPREAD {
            return Ok(None);
        }

        let pos0 = median_u64(&mut a_pos);
        let mate0 = median_u64(&mut b_pos);

        return Ok(Some(SvCall {
            chrom: chrom.into(),
            pos: pos0 + 1,
            end: pos0 + 1,
            svtype: "BND".into(),
            len: 0,
            score: a_pos.len() as i32,
            reason: "PE".into(),
            mate_chrom: mate_chrom_name,
            mate_pos: Some(mate0 + 1),
        }));
    }

    // -------- INS from SR-like soft clips --------
    if hs.reason == "soft_clip_sr" {
        const MIN_MAPQ: u8 = 20;
        const MIN_SR_SUPPORT: usize = 4;
        const MAX_POS_SPREAD: u64 = 150;
        const MIN_INS_LEN: isize = 50;

        let min_clip_sr: isize = (k as isize).max(20);

        let mut clip_lens: Vec<isize> = Vec::new();
        let mut clip_pos0: Vec<u64> = Vec::new();
        let mut seen_qname: HashSet<Vec<u8>> = HashSet::new();

        for r in &support {
            if r.is_secondary() || r.is_duplicate() {
                continue;
            }
            if r.mapq() < MIN_MAPQ {
                continue;
            }

            let split_like = r.is_supplementary() || r.aux(b"SA").is_ok();
            if !split_like {
                continue;
            }

            let qn = r.qname().to_vec();
            if !seen_qname.insert(qn) {
                continue;
            }

            let pos0 = r.pos().max(0) as u64;
            let cig = r.cigar();

            // check softclip on either end
            if let Some(Cigar::SoftClip(n)) = cig.iter().next() {
                let n = *n as isize;
                if n >= min_clip_sr {
                    clip_lens.push(n);
                    clip_pos0.push(pos0);
                    continue;
                }
            }
            if let Some(Cigar::SoftClip(n)) = cig.iter().last() {
                let n = *n as isize;
                if n >= min_clip_sr {
                    clip_lens.push(n);
                    clip_pos0.push(pos0 + ref_consumed_len(&r.cigar()));
                    continue;
                }
            }
        }

        if clip_lens.len() < MIN_SR_SUPPORT {
            return Ok(None);
        }

        clip_pos0.sort_unstable();
        let spread = clip_pos0[clip_pos0.len() - 1] - clip_pos0[0];
        if spread > MAX_POS_SPREAD {
            return Ok(None);
        }

        let approx_len = {
            clip_lens.sort_unstable();
            clip_lens[clip_lens.len() / 2]
        };
        if approx_len < MIN_INS_LEN {
            return Ok(None);
        }

        let ins0 = median_u64(&mut clip_pos0);
        let pos_vcf = ins0 + 1;

        return Ok(Some(SvCall {
            chrom: chrom.into(),
            pos: pos_vcf,
            end: pos_vcf,
            svtype: "INS".into(),
            len: approx_len,
            score: clip_lens.len() as i32,
            reason: "SR".into(),
            mate_chrom: None,
            mate_pos: None,
        }));
    }

    Ok(None)
}


