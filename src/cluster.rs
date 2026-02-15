use crate::confirm::SvCall;

pub fn cluster_calls(mut calls: Vec<SvCall>, radius: u64) -> Vec<SvCall> {
    calls.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then_with(|| a.svtype.cmp(&b.svtype))
            .then_with(|| a.pos.cmp(&b.pos))
            .then_with(|| a.end.cmp(&b.end))
    });

    let mut out: Vec<SvCall> = Vec::new();
    let mut cur: Option<SvCall> = None;

    for c in calls.into_iter() {
        if let Some(ref mut x) = cur {
            let same = x.chrom == c.chrom && x.svtype == c.svtype;
            let pos_close = x.pos.abs_diff(c.pos) <= radius;
            let end_close = x.end.abs_diff(c.end) <= radius;

            if same && pos_close && end_close {
                // Lumpy-ish merge behaviour: expand interval rather than averaging
                x.pos = x.pos.min(c.pos);
                x.end = x.end.max(c.end);

                // support score: keep max
                x.score = x.score.max(c.score);

                // length: recompute for DEL if possible
                if x.svtype == "DEL" {
                    let len = x.end.saturating_sub(x.pos) as isize;
                    x.len = -(len as isize);
                }

                // merge reasons
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
