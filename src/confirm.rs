use anyhow::Result;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;
use serde::{Deserialize, Serialize};

use crate::hotspot::Hotspot;
use crate::suffix::SuffixWorkspace;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SvCall {
    pub chrom: String,
    pub pos: u64,     // VCF POS (1-based when written)
    pub end: u64,     // VCF END (1-based conceptually; we store as 1-based coordinate)
    pub svtype: String,
    pub len: isize,
    pub score: i32,
    pub reason: String,
}

// median for u64
fn median_u64(v: &mut Vec<u64>) -> u64 {
    v.sort_unstable();
    v[v.len() / 2]
}

// median for isize
fn median_isize(v: &mut Vec<isize>) -> isize {
    v.sort_unstable();
    v[v.len() / 2]
}

pub fn confirm_breakpoint(
    bin_start: u64,
    chrom: &str,
    ref_window: &[u8],
    _sfx: &SuffixWorkspace,
    hs: &Hotspot,
    band: usize,
    insert_mean: f64,
    reads: &[bam::Record],
) -> Result<Option<SvCall>> {
    let center = hs.local_pos.min(ref_window.len().saturating_sub(1));
    let local_start = center.saturating_sub(band);
    let local_end = (center + band).min(ref_window.len());

    // NOTE: BAM positions are 0-based. VCF POS is 1-based.
    let g_start0 = bin_start + local_start as u64; // 0-based
    let g_end0 = bin_start + local_end as u64;     // 0-based end

    // Collect reads overlapping local window
    let mut support: Vec<&bam::Record> = Vec::new();
    for r in reads.iter() {
        let pos0 = r.pos().max(0) as u64;
        let end0 = pos0 + r.seq_len() as u64;
        if end0 > g_start0 && pos0 < g_end0 {
            support.push(r);
        }
    }
    if support.is_empty() {
        return Ok(None);
    }

    // -------------------------
    // DEL: discordant pairs -> estimate start/end from anchors
    // -------------------------
    if hs.reason == "discordant" {
        const MIN_TLEN_DEV: f64 = 200.0; // how much larger than mean
        const MIN_DEL_BP: u64 = 50;      // don't emit micro events

        let mut left_bps: Vec<u64> = Vec::new();
        let mut right_bps: Vec<u64> = Vec::new();

        for r in &support {
            if !(r.is_paired() && !r.is_unmapped() && !r.is_mate_unmapped()) {
                continue;
            }
            if r.tid() != r.mtid() {
                continue; // not handling BND here (yet)
            }

            let tlen = r.insert_size().abs() as f64;
            if tlen - insert_mean <= MIN_TLEN_DEV {
                continue;
            }

            // Anchor positions (0-based)
            let rpos0 = r.pos().max(0) as u64;
            let mpos0 = r.mpos().max(0) as u64;

            let left0 = rpos0.min(mpos0);
            let right0 = rpos0.max(mpos0);

            // For deletions, breakpoint tends to be near end of left read and start of right read.
            // Left breakpoint approx = left + read_len
            let read_len = r.seq_len() as u64;
            let left_bp0 = left0 + read_len;
            let right_bp0 = right0;

            if right_bp0 > left_bp0 {
                left_bps.push(left_bp0);
                right_bps.push(right_bp0);
            }
        }

        if left_bps.is_empty() || right_bps.is_empty() {
            return Ok(None);
        }

        let start0 = median_u64(&mut left_bps);
        let end0 = median_u64(&mut right_bps);

        if end0 <= start0 || (end0 - start0) < MIN_DEL_BP {
            return Ok(None);
        }

        let svlen = -((end0 - start0) as isize);

        // Convert to VCF-style coordinates:
        // - VCF POS is 1-based, so POS = start0 + 1
        // - END is typically 1-based inclusive-ish in many callers; Lumpy uses END as coordinate.
        // We'll store END as end0 (0-based) + 1.
        let pos_vcf = start0 + 1;
        let end_vcf = end0 + 1;

        return Ok(Some(SvCall {
            chrom: chrom.into(),
            pos: pos_vcf,
            end: end_vcf,
            svtype: "DEL".into(),
            len: svlen,
            score: (left_bps.len() as i32),
            reason: hs.reason.clone(),
        }));
    }

    // -------------------------
    // INS: soft clips -> estimate inserted length from clip size (approx)
    // -------------------------
    if hs.reason == "soft_clip" {
        const MIN_SC_SUPPORT: usize = 2;
        const MIN_CLIP: isize = 15;

        let mut clip_lens: Vec<isize> = Vec::new();
        let mut clip_positions0: Vec<u64> = Vec::new();

        for r in &support {
            let pos0 = r.pos().max(0) as u64;
            let cig = r.cigar();

            if let Some(Cigar::SoftClip(n)) = cig.iter().next() {
                clip_lens.push(*n as isize);
                clip_positions0.push(pos0);
            } else if let Some(Cigar::SoftClip(n)) = cig.iter().last() {
                clip_lens.push(*n as isize);
                clip_positions0.push(pos0 + r.seq_len() as u64);
            }
        }

        if clip_lens.len() < MIN_SC_SUPPORT {
            return Ok(None);
        }

        let approx_len = median_isize(&mut clip_lens);
        if approx_len < MIN_CLIP {
            return Ok(None);
        }

        // choose a representative position near the clip cluster
        let ins_pos0 = median_u64(&mut clip_positions0);

        let pos_vcf = ins_pos0 + 1;
        let end_vcf = pos_vcf; // INS: END often equals POS unless assembled

        return Ok(Some(SvCall {
            chrom: chrom.into(),
            pos: pos_vcf,
            end: end_vcf,
            svtype: "INS".into(),
            len: approx_len,
            score: clip_lens.len() as i32,
            reason: hs.reason.clone(),
        }));
    }

    Ok(None)
}
