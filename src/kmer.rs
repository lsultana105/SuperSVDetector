use hashbrown::HashMap;
use rust_htslib::bam;
use rust_htslib::bam::record::{Cigar, Seq};

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

pub fn discover_hotspots(
    reads: &[bam::Record],
    mindex: &MinimizerIndex,
    k: usize,
    bin_start: u64,
    insert_mean: f64,
    insert_sd: f64,
) -> Vec<Hotspot> {
    let mut evidence: HashMap<usize, EvidenceCounts> = HashMap::new();

    // ---------- tuning knobs (aiming closer to Lumpy) ----------
    const MIN_MAPQ: u8 = 20;

    // Lumpy-ish discordant filter (discordant_z≈5)
    const Z_CUTOFF: f64 = 5.0;

    // Hotspot support thresholds (raise/lower as needed)
    const MIN_DEL_PE_SUPPORT: usize = 6;
    const MIN_BND_SUPPORT: usize = 3;

    // IMPORTANT: SR evidence must be stronger than “any soft clip”
    const MIN_SR_SUPPORT: usize = 4;

    // window sizes
    const WIN_DISCORDANT: usize = 500;
    const WIN_BND: usize = 500;
    const WIN_SR: usize = 100;
    const WIN_RARE: usize = 200;

    // Rare kmers disabled
    const MIN_RARE_SUPPORT: usize = 999999;

    let min_clip = (k as u32).max(20);

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
        let local = (gpos0 - bin_start) as usize;

        // ---------------------------
        // 1) PE evidence
        // ---------------------------
        if r.is_paired() && !r.is_unmapped() && !r.is_mate_unmapped() {
            // inter-chrom = BND-ish signal
            if r.tid() != r.mtid() {
                // count one per pair
                if r.pos() < r.mpos() {
                    let b = (local / WIN_BND) * WIN_BND;
                    evidence.entry(b).or_default().bnd_pe += 1;
                }
                continue;
            }

            // same chrom: discordant DEL-like
            if !r.is_proper_pair() && insert_sd > 0.0 && r.pos() < r.mpos() {
                let tlen = r.insert_size().abs() as f64;
                let z = (tlen - insert_mean) / insert_sd;
                if z >= Z_CUTOFF {
                    let b = (local / WIN_DISCORDANT) * WIN_DISCORDANT;
                    evidence.entry(b).or_default().discordant_del += 1;
                    continue;
                }
            }
        }

        // ---------------------------
        // 2) SR-like evidence (soft-clip but only if split-ish)
        //    Lumpy uses splitters BAM. Approximate that by requiring:
        //      - supplementary OR SA tag present
        //      - and a sufficiently long soft clip
        // ---------------------------
        let has_sa = r.aux(b"SA").is_ok();
        let is_split_like = has_sa || r.is_supplementary();

        if is_split_like {
            let cig = r.cigar();
            let mut ok = false;

            if let Some(Cigar::SoftClip(n)) = cig.iter().next() {
                if *n >= min_clip { ok = true; }
            }
            if !ok {
                if let Some(Cigar::SoftClip(n)) = cig.iter().last() {
                    if *n >= min_clip { ok = true; }
                }
            }

            if ok {
                let b = (local / WIN_SR) * WIN_SR;
                evidence.entry(b).or_default().soft_clip_sr += 1;
                continue;
            }
        }

        // ---------------------------
        // 3) Rare kmers (disabled)
        // ---------------------------
        if r.seq_len() >= k {
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
                let b = (local / WIN_RARE) * WIN_RARE;
                evidence.entry(b).or_default().rare_kmer += 1;
            }
        }
    }

    let mut out: Vec<Hotspot> = Vec::new();
    for (bin, c) in evidence.into_iter() {
        if c.discordant_del >= MIN_DEL_PE_SUPPORT {
            out.push(Hotspot { local_pos: bin, reason: "discordant_del".into() });
        }
        if c.bnd_pe >= MIN_BND_SUPPORT {
            out.push(Hotspot { local_pos: bin, reason: "bnd_pe".into() });
        }
        if c.soft_clip_sr >= MIN_SR_SUPPORT {
            out.push(Hotspot { local_pos: bin, reason: "soft_clip_sr".into() });
        }
        if c.rare_kmer >= MIN_RARE_SUPPORT {
            out.push(Hotspot { local_pos: bin, reason: "rare_kmer".into() });
        }
    }

    out.sort_by_key(|h| h.local_pos);
    out
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