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
    let mut current_clustered_call: Option<SvCall> = None;

    for call in calls.into_iter() {
        if let Some(ref mut current_call) = current_clustered_call {
            let same = current_call.chrom == call.chrom && current_call.svtype == call.svtype;
            let pos_close = current_call.pos.abs_diff(call.pos) <= radius;
            let end_close = current_call.end.abs_diff(call.end) <= radius;

            if same && pos_close && end_close {
                // Lumpy-ish merge behaviour: expand interval rather than averaging
                current_call.pos = current_call.pos.min(call.pos);
                current_call.end = current_call.end.max(call.end);

                // support score: keep max
                current_call.score = current_call.score.max(call.score);

                // length: recompute for DEL if possible
                if current_call.svtype == "DEL" {
                    let len = current_call.end.saturating_sub(current_call.pos) as isize;
                    current_call.len = -(len as isize);
                }

                // merge reasons
                if !current_call.reason.contains(&call.reason) {
                    if !current_call.reason.is_empty() {
                        current_call.reason.push(',');
                    }
                    current_call.reason.push_str(&call.reason);
                }
                continue;
            }

            out.push(current_clustered_call.take().unwrap());
        }
        current_clustered_call = Some(call);
    }

    if let Some(x) = current_clustered_call {
        out.push(x);
    }
    out
}
