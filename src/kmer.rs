use std::collections::{HashMap, HashSet};

use rust_htslib::bam;
use rust_htslib::bam::record::Seq;

use crate::hotspot::Hotspot;

pub struct MinimizerIndex {
    pub map: HashMap<u64, Vec<usize>>,
}

impl MinimizerIndex {
    pub fn build(ref_window: &[u8], k: usize, w: usize) -> Self {
        assert!(w >= k);
        let mut map: HashMap<u64, Vec<usize>> = HashMap::new();
        if ref_window.len() < w {
            return Self { map };
        }

        let enc = |b: u8| match b {
            b'A' | b'a' => Some(0u64),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' => Some(3),
            _ => None,
        };

        for i in 0..=ref_window.len() - w {
            let mut best: Option<(u64, usize)> = None;

            for j in i..=i + w - k {
                let mut val = 0u64;
                let mut ok = true;

                for t in 0..k {
                    if let Some(x) = enc(ref_window[j + t]) {
                        val = (val << 2) | x;
                    } else {
                        ok = false;
                        break;
                    }
                }

                if ok && best.map_or(true, |(v, _)| val < v) {
                    best = Some((val, j));
                }
            }

            if let Some((val, jpos)) = best {
                map.entry(val).or_default().push(jpos);
            }
        }

        Self { map }
    }
}

#[derive(Default, Clone, Debug)]
struct EvidenceCounts {
    discordant_del: usize,
    bnd_pe: usize,
    soft_clip_sr: usize,
    rare_kmer: usize,
}

fn has_sa_tag(r: &bam::Record) -> bool {
    r.aux(b"SA").is_ok()
}

fn collapse_hotspots(mut raw: Vec<(usize, String, usize)>) -> Vec<Hotspot> {
    if raw.is_empty() {
        return Vec::new();
    }

    raw.sort_by(|a, b| {
        let ord = a.1.cmp(&b.1);
        if ord == std::cmp::Ordering::Equal {
            a.0.cmp(&b.0)
        } else {
            ord
        }
    });

    let mut out: Vec<Hotspot> = Vec::new();

    let mut cur_pos = raw[0].0;
    let mut cur_reason = raw[0].1.clone();
    let mut cur_support = raw[0].2;

    const MERGE_WIN: usize = 300;

    for (pos, reason, support) in raw.into_iter().skip(1) {
        if reason == cur_reason && pos.abs_diff(cur_pos) <= MERGE_WIN {
            if support > cur_support {
                cur_pos = pos;
                cur_support = support;
            }
        } else {
            out.push(Hotspot {
                local_pos: cur_pos,
                reason: cur_reason.clone(),
            });
            cur_pos = pos;
            cur_reason = reason;
            cur_support = support;
        }
    }

    out.push(Hotspot {
        local_pos: cur_pos,
        reason: cur_reason,
    });

    out.sort_by_key(|h| h.local_pos);
    out
}

pub fn discover_hotspots(
    reads: &[bam::Record],
    mindex: &MinimizerIndex,
    k: usize,
    bin_start: u64,
    insert_mean: f64,
    insert_sd: f64,
) -> Vec<Hotspot> {
    use rust_htslib::bam::record::Cigar;

    let mut evidence: HashMap<usize, EvidenceCounts> = HashMap::new();

    const MIN_MAPQ: u8 = 20;
    const Z_CUTOFF: f64 = 3.0;

    // Make deletion hotspot generation permissive again.
    // Let confirm.rs be the stricter stage.
    const MIN_DEL_PE_SUPPORT: usize = 2;
    const MIN_BND_SUPPORT: usize = 3;
    const MIN_SR_SUPPORT: usize = 3;

    const WIN_DISCORDANT: usize = 500;
    const WIN_BND: usize = 500;
    const WIN_SR: usize = 100;
    const WIN_RARE: usize = 200;

    const MIN_RARE_SUPPORT: usize = 999999;

    let min_clip_sr: u32 = 20;

    let hard_tlen_min: f64 = if insert_sd > 0.0 {
        (insert_mean + 3.0 * insert_sd).max(800.0)
    } else {
        800.0
    };

    let mut seen_del_pe: HashSet<Vec<u8>> = HashSet::new();
    let mut seen_bnd_pe: HashSet<Vec<u8>> = HashSet::new();
    let mut seen_sr: HashSet<Vec<u8>> = HashSet::new();

    for r in reads {
        if r.is_secondary() || r.is_duplicate() {
            continue;
        }
        if r.mapq() < MIN_MAPQ {
            continue;
        }

        let gpos0 = r.pos().max(0) as u64;
        if gpos0 < bin_start {
            continue;
        }

        if r.is_paired()
            && !r.is_unmapped()
            && !r.is_mate_unmapped()
            && !r.is_supplementary()
        {
            // Inter-chromosomal pairs -> BND-like evidence
            if r.tid() != r.mtid() {
                let qn = r.qname().to_vec();
                if !seen_bnd_pe.insert(qn) {
                    continue;
                }

                let local = (gpos0 - bin_start) as usize;
                let b = (local / WIN_BND) * WIN_BND;
                evidence.entry(b).or_default().bnd_pe += 1;
                continue;
            }

            let tlen = r.insert_size().abs() as f64;
            let mut is_del_like = false;

            if insert_sd > 0.0 {
                let z = (tlen - insert_mean) / insert_sd;
                if z >= Z_CUTOFF {
                    is_del_like = true;
                }
            }

            if tlen >= hard_tlen_min {
                is_del_like = true;
            }

            if is_del_like {
                let qn = r.qname().to_vec();
                if !seen_del_pe.insert(qn) {
                    continue;
                }

                let rpos0 = r.pos().max(0) as u64;
                let mpos0 = r.mpos().max(0) as u64;
                let left0 = rpos0.min(mpos0);

                if left0 >= bin_start {
                    let left_local = (left0 - bin_start) as usize;
                    let b = (left_local / WIN_DISCORDANT) * WIN_DISCORDANT;
                    evidence.entry(b).or_default().discordant_del += 1;
                }
                continue;
            }
        }

        let qn = r.qname().to_vec();
        if seen_sr.contains(&qn) {
            continue;
        }

        // Soft-clipped split-read-like evidence
        let cig = r.cigar();
        let mut clip_ok = false;

        if let Some(Cigar::SoftClip(n)) = cig.iter().next() {
            if *n >= min_clip_sr {
                clip_ok = true;
            }
        }
        if !clip_ok {
            if let Some(Cigar::SoftClip(n)) = cig.iter().last() {
                if *n >= min_clip_sr {
                    clip_ok = true;
                }
            }
        }

        // Keep the stricter SA/supplementary requirement for SR evidence
        if clip_ok && (has_sa_tag(r) || r.is_supplementary()) {
            seen_sr.insert(qn);
            let local = (gpos0 - bin_start) as usize;
            let b = (local / WIN_SR) * WIN_SR;
            evidence.entry(b).or_default().soft_clip_sr += 1;
            continue;
        }

        // Rare k-mer evidence (disabled by default)
        if MIN_RARE_SUPPORT != 999999 && r.seq_len() >= k {
            let seq = r.seq();
            let (mut rare, mut total) = (0usize, 0usize);

            for i in (0..=seq.len() - k).step_by(5) {
                if let Some(v) = encode_kmer(&seq, i, k) {
                    total += 1;
                    if mindex.map.get(&v).is_none() {
                        rare += 1;
                    }
                }
            }

            if total >= 5 && rare * 2 > total {
                let local = (gpos0 - bin_start) as usize;
                let b = (local / WIN_RARE) * WIN_RARE;
                evidence.entry(b).or_default().rare_kmer += 1;
            }
        }
    }

    let mut raw: Vec<(usize, String, usize)> = Vec::new();

    for (bin, c) in evidence.into_iter() {
        if c.discordant_del >= MIN_DEL_PE_SUPPORT {
            raw.push((bin, "discordant_del".into(), c.discordant_del));
        }
        if c.bnd_pe >= MIN_BND_SUPPORT {
            raw.push((bin, "bnd_pe".into(), c.bnd_pe));
        }
        if c.soft_clip_sr >= MIN_SR_SUPPORT {
            raw.push((bin, "soft_clip_sr".into(), c.soft_clip_sr));
        }
        if c.rare_kmer >= MIN_RARE_SUPPORT {
            raw.push((bin, "rare_kmer".into(), c.rare_kmer));
        }
    }

    collapse_hotspots(raw)
}

fn encode_kmer(seq: &Seq, i: usize, k: usize) -> Option<u64> {
    let mut val = 0u64;
    let bytes = seq.as_bytes();
    for t in 0..k {
        val = (val << 2)
            | match bytes[i + t] {
            1 | b'A' | b'a' => 0,
            2 | b'C' | b'c' => 1,
            4 | b'G' | b'g' => 2,
            8 | b'T' | b't' => 3,
            _ => return None,
        };
    }
    Some(val)
}