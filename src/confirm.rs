use anyhow::Result;
use bio::alignment::pairwise::{Aligner, Scoring};
use bio::alignment::AlignmentOperation;
use serde::{Serialize, Deserialize};

use crate::bins::Bin;
use crate::kmer::Hotspot;
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

pub fn confirm_breakpoint(
    bin: &Bin,
    ref_window: &[u8],
    sfx: &SuffixWorkspace,
    hs: &Hotspot,
    mut band: usize,
) -> Result<Option<SvCall>> {
    let center = hs.local_pos.min(ref_window.len().saturating_sub(1));

    // --- Adjust band size depending on hotspot reason ---
    if hs.reason.contains("discordant_large_insert") {
        band = band.max(128); // allow bigger gaps for deletions
    } else if hs.reason.contains("discordant_small_insert") {
        band = band.max(64); // allow tighter gaps for small indels
    } else if hs.reason.contains("discordant_orientation") {
        band = band.max(96); // orientation anomalies may imply inversions
    }

    // --- Banded alignment ---
    let left = center.saturating_sub(band);
    let right = (center + band).min(ref_window.len());
    if right <= left + 20 {
        return Ok(None);
    }
    let query = &ref_window[left..right];

    let scoring = Scoring::new(-6, -1, |a: u8, b: u8| if a == b { 2 } else { -2 });
    let mut aligner = Aligner::with_capacity_and_scoring(ref_window.len(), query.len(), scoring);
    let aln = aligner.global(ref_window, query);

    // --- Count longest insertion/deletion runs ---
    let mut ins = 0isize;
    let mut del = 0isize;
    let mut cur_ins = 0isize;
    let mut cur_del = 0isize;
    for op in aln.operations.iter() {
        match op {
            AlignmentOperation::Ins => {
                cur_ins += 1;
                ins = ins.max(cur_ins);
                cur_del = 0;
            }
            AlignmentOperation::Del => {
                cur_del += 1;
                del = del.max(cur_del);
                cur_ins = 0;
            }
            _ => {
                cur_ins = 0;
                cur_del = 0;
            }
        }
    }

    // --- Decide SV type ---
    let (svtype, len) = if del >= 10 {
        ("DEL".to_string(), -(del as isize))
    } else if ins >= 10 {
        ("INS".to_string(), (ins as isize))
    } else {
        return Ok(None);
    };

    let pos = bin.start + center as u64 + 1;
    Ok(Some(SvCall {
        chrom: bin.chrom.clone(),
        pos,
        svtype,
        len,
        score: aln.score,
        reason: hs.reason.clone(), // keep hotspot reason for traceability
    }))
}
