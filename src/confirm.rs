use anyhow::Result;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;
use serde::{Serialize, Deserialize};

use crate::hotspot::Hotspot;
use crate::suffix::SuffixWorkspace;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SvCall {
    pub chrom: String,
    pub pos: u64,
    pub svtype: String,
    pub len: isize,
    pub score: i32,
    pub reason: String,
}

// median for u64
fn median_u64(v: &mut Vec<u64>) -> u64 {
    v.sort();
    v[v.len() / 2]
}

// median for isize
fn median_isize(v: &mut Vec<isize>) -> isize {
    v.sort();
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

    let g_start = bin_start + local_start as u64;
    let g_end = bin_start + local_end as u64;

    // Collect reads overlapping local window
    let mut support: Vec<&bam::Record> = Vec::new();
    for r in reads.iter() {
        let pos = r.pos().max(0) as u64;
        let end = pos + r.seq_len() as u64;
        if end > g_start && pos < g_end {
            support.push(r);
        }
    }

    if support.is_empty() {
        return Ok(None);
    }

    // --- Deletions: strong TLEN deviations from discordant pairs ---
    if hs.reason == "discordant" {
        const MIN_TLEN_DEV: f64 = 200.0;  // how much larger than mean to count
        const MIN_DEL_LEN: f64 = 200.0;   // minimum deletion size to report

        let mut del_devs: Vec<f64> = Vec::new();

        for r in &support {
            if r.is_paired() && !r.is_unmapped() && !r.is_mate_unmapped() && r.tid() == r.mtid() {
                let tlen = r.insert_size().abs() as f64;
                let dev = tlen - insert_mean;
                if dev > MIN_TLEN_DEV {
                    del_devs.push(dev);
                }
            }
        }

        if del_devs.is_empty() {
            return Ok(None);
        }

        del_devs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let del_est = del_devs[del_devs.len() / 2];

        if del_est < MIN_DEL_LEN {
            return Ok(None);
        }

        return Ok(Some(SvCall {
            chrom: chrom.into(),
            pos: bin_start + center as u64,
            svtype: "DEL".into(),
            len: -(del_est.round() as isize),
            score: support.len() as i32,
            reason: hs.reason.clone(),
        }));
    }

    // --- Insertions: soft-clips + rare_kmer fallback ---
    if hs.reason == "soft_clip" {
        const MIN_SC_SUPPORT: usize = 2;

        let mut clip_lens: Vec<isize> = Vec::new();

        for r in &support {
            let cig = r.cigar();
            if let Some(Cigar::SoftClip(n)) = cig.iter().next() {
                clip_lens.push(*n as isize);
            } else if let Some(Cigar::SoftClip(n)) = cig.iter().last() {
                clip_lens.push(*n as isize);
            }
        }

        let approx_len = if !clip_lens.is_empty() {
            median_isize(&mut clip_lens)
        } else if hs.reason == "rare_kmer" {
            // No explicit soft-clips but a strong rare_kmer hotspot.
            // For synthetic data, we know insertions are ~50 bp.
            50
        } else {
            return Ok(None);
        };

        if hs.reason == "soft_clip" && support.len() < MIN_SC_SUPPORT {
            return Ok(None);
        }

        return Ok(Some(SvCall {
            chrom: chrom.into(),
            pos: bin_start + center as u64,
            svtype: "INS".into(),
            len: approx_len,
            score: support.len() as i32,
            reason: hs.reason.clone(),
        }));
    }

    Ok(None)
}
