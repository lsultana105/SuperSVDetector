//! Simple spatial clustering of SV calls across bins.
//! Groups calls by (chrom, svtype) and merges if positions are within `radius`.

use crate::confirm::SvCall;

pub fn cluster_calls(mut calls: Vec<SvCall>, radius: u64) -> Vec<SvCall> {
    // Sort calls by (chrom, svtype, pos)
    calls.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then_with(|| a.svtype.cmp(&b.svtype))
            .then_with(|| a.pos.cmp(&b.pos))
    });

    let mut out: Vec<SvCall> = Vec::new();
    let mut cur: Option<SvCall> = None;

    for c in calls.into_iter() {
        if let Some(ref mut x) = cur {
            if x.chrom == c.chrom
                && x.svtype == c.svtype
                && c.pos.saturating_sub(x.pos) <= radius
            {
                // --- merge logic ---
                // average position/length
                x.pos = (x.pos + c.pos) / 2;
                x.len = ((x.len + c.len) / 2) as isize;
                // keep max score
                x.score = x.score.max(c.score);
                // merge reasons (avoid duplicates)
                if !x.reason.contains(&c.reason) {
                    if !x.reason.is_empty() {
                        x.reason.push(',');
                    }
                    x.reason.push_str(&c.reason);
                }
                continue;
            }
            out.push(cur.take().unwrap());
        }
        cur = Some(c);
    }

    if let Some(x) = cur {
        out.push(x);
    }
    out
}
