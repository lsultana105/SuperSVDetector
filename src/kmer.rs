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

pub fn discover_hotspots(reads: &[bam::Record], mindex: &MinimizerIndex, k: usize) -> Vec<Hotspot> {
    let mut evidence_map: HashMap<usize, (usize, String)> = HashMap::new();

    const MIN_RP_SUPPORT: usize = 4;
    const MIN_SC_SUPPORT: usize = 1;    // allow even 1 soft-clipped read
    const MIN_RARE_SUPPORT: usize = 9999; // effectively disable rare_kmer hotspots for now

    println!("--- Starting Hotspot Discovery ---");
    println!("Total reads to process: {}", reads.len());

    for r in reads {
        if r.is_secondary() || r.is_duplicate() {
            continue;
        }

        let pos = r.pos().max(0) as usize;
        let mut evidence_found = false;
        let mut reason = "";

        // 1. Discordant Pairs → candidate deletions
        if r.is_paired() && !r.is_unmapped() && !r.is_mate_unmapped() {
            let isize = r.insert_size().abs() as usize;
            if isize > 1000 || isize < 100 || r.is_reverse() == r.is_mate_reverse() {
                evidence_found = true;
                reason = "discordant";
            }
        }

        // 2. Soft-clips → strong insertion / breakpoint evidence
        if !evidence_found {
            let cigar = r.cigar();
            if let Some(Cigar::SoftClip(n)) = cigar.iter().next() {
                if *n as usize >= k {
                    evidence_found = true;
                    reason = "soft_clip";
                }
            } else if let Some(Cigar::SoftClip(n)) = cigar.iter().last() {
                if *n as usize >= k {
                    evidence_found = true;
                    reason = "soft_clip";
                }
            }
        }

        // 3. Rare Minimizers → weird / novel sequence (fallback signal)
        if !evidence_found && r.seq_len() >= k {
            let seq = r.seq();
            let (mut rare, mut total) = (0, 0);
            for i in (0..=seq.len() - k).step_by(5) {
                if let Some(v) = encode_kmer(&seq, i, k) {
                    total += 1;
                    if mindex.map.get(&v).is_none() {
                        rare += 1;
                    }
                }
            }
            if total >= 5 && rare * 2 > total {
                evidence_found = true;
                reason = "rare_kmer";
            }
        }

        if evidence_found {
            // Group by 100bp window to count supporting reads
            let bin = (pos / 100) * 100;
            let entry = evidence_map.entry(bin).or_insert((0, reason.to_string()));
            entry.0 += 1;
        }
    }

    println!("Candidate sites found before filtering: {}", evidence_map.len());

    let mut final_hs: Vec<Hotspot> = evidence_map
        .into_iter()
        .filter(|(_, (count, reason))| match reason.as_str() {
            "soft_clip"   => *count >= MIN_SC_SUPPORT,
            "discordant"  => *count >= MIN_RP_SUPPORT,
            "rare_kmer"   => *count >= MIN_RARE_SUPPORT,
            _ => false,
        })
        .map(|(bin, (_, reason))| Hotspot {
            local_pos: bin,
            reason,
        })
        .collect();

    println!("Final High-Confidence Hotspots: {}", final_hs.len());

    final_hs.sort_by_key(|h| h.local_pos);
    final_hs
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
