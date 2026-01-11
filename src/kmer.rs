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
        if ref_window.len() < w { return Self { map }; }

        let enc = |b: u8| match b {
            b'A' | b'a' => Some(0u64), b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2), b'T' | b't' => Some(3), _ => None,
        };

        for i in 0..=ref_window.len() - w {
            let mut best: Option<(u64, usize)> = None;
            for j in i..=i + w - k {
                let mut val = 0u64;
                let mut ok = true;
                for t in 0..k {
                    if let Some(x) = enc(ref_window[j + t]) { val = (val << 2) | x; }
                    else { ok = false; break; }
                }
                if ok && best.map_or(true, |(v, _)| val < v) { best = Some((val, j)); }
            }
            if let Some((val, jpos)) = best { map.entry(val).or_default().push(jpos); }
        }
        Self { map }
    }
}

pub fn discover_hotspots(reads: &[bam::Record], mindex: &MinimizerIndex, k: usize) -> Vec<Hotspot> {
    let mut hs = Vec::new();
    for r in reads {
        let pos = r.pos().max(0) as usize;
        let mut flagged = false;

        // 1. Quick Check: Discordant Pairs
        if r.is_paired() && !r.is_unmapped() && !r.is_mate_unmapped() {
            let isize = r.insert_size().abs() as usize;
            if isize > 1000 || isize < 100 || r.is_reverse() == r.is_mate_reverse() {
                hs.push(Hotspot { local_pos: pos, reason: "discordant".into() });
                flagged = true;
            }
        }

        // 2. Medium Check: Soft-clips
        if !flagged {
            let cigar = r.cigar();
            if let Some(Cigar::SoftClip(n)) = cigar.iter().next() {
                if *n as usize >= k { hs.push(Hotspot { local_pos: pos, reason: "soft_L".into() }); flagged = true; }
            }
            if !flagged {
                if let Some(Cigar::SoftClip(n)) = cigar.iter().last() {
                    if *n as usize >= k { hs.push(Hotspot { local_pos: pos + r.seq_len(), reason: "soft_R".into() }); flagged = true; }
                }
            }
        }

        // 3. Heavy Check: Minimizer Rarity (Only if no other clues found)
        if !flagged && r.seq_len() >= k {
            let seq = r.seq();
            let (mut rare, mut total) = (0, 0);
            for i in 0..=seq.len() - k {
                if let Some(v) = encode_kmer(&seq, i, k) {
                    total += 1;
                    if mindex.map.get(&v).map_or(0, |p| p.len()) == 0 { rare += 1; }
                }
            }
            if total >= 10 && rare * 2 > total { hs.push(Hotspot { local_pos: pos, reason: "rare".into() }); }
        }
    }
    hs.sort_by_key(|h| h.local_pos);
    hs.dedup_by(|a, b| (a.local_pos / 50) == (b.local_pos / 50));
    hs
}

fn encode_kmer(seq: &Seq, i: usize, k: usize) -> Option<u64> {
    let mut val = 0u64;
    let bytes = seq.as_bytes();
    for t in 0..k {
        val = (val << 2) | match bytes[i + t] {
            1 | b'A' | b'a' => 0, 2 | b'C' | b'c' => 1,
            4 | b'G' | b'g' => 2, 8 | b'T' | b't' => 3, _ => return None,
        };
    }
    Some(val)
}