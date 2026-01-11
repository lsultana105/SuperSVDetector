use anyhow::Result;
use bio::alignment::pairwise::{Aligner, Scoring};
use bio::alignment::AlignmentOperation;
use crate::suffix::SuffixWorkspace;
use crate::hotspot::Hotspot;
use serde::{Serialize, Deserialize};

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
    bin_start: u64,
    chrom: &str,
    ref_window: &[u8],
    _sfx: &SuffixWorkspace,
    hs: &Hotspot,
    band: usize,
) -> Result<Option<SvCall>> {
    let center = hs.local_pos.min(ref_window.len().saturating_sub(1));

    // REDUCE SEARCH SPACE: Only align a small local patch
    let search_radius = band.max(200);
    let ref_start = center.saturating_sub(search_radius);
    let ref_end = (center + search_radius).min(ref_window.len());

    if ref_end <= ref_start + 20 { return Ok(None); }

    let local_ref = &ref_window[ref_start..ref_end];
    let query_size = 100.min(local_ref.len() / 2); // Small probe
    let query = &local_ref[local_ref.len()/2 - query_size/2 .. local_ref.len()/2 + query_size/2];

    let scoring = Scoring::new(-6, -1, |a: u8, b: u8| if a == b { 2 } else { -2 });
    let mut aligner = Aligner::with_capacity_and_scoring(local_ref.len(), query.len(), scoring);
    let aln = aligner.global(local_ref, query);

    let (mut ins, mut del, mut c_ins, mut c_del) = (0, 0, 0, 0);
    for op in aln.operations.iter() {
        match op {
            AlignmentOperation::Ins => { c_ins += 1; ins = ins.max(c_ins); c_del = 0; }
            AlignmentOperation::Del => { c_del += 1; del = del.max(c_del); c_ins = 0; }
            _ => { c_ins = 0; c_del = 0; }
        }
    }

    let (svtype, len) = if del >= 10 { ("DEL".into(), -(del as isize)) }
    else if ins >= 10 { ("INS".into(), ins as isize) }
    else { return Ok(None); };

    Ok(Some(SvCall {
        chrom: chrom.into(),
        pos: bin_start + center as u64,
        svtype, len, score: aln.score, reason: hs.reason.clone(),
    }))
}