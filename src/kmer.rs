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

/// Evidence counts per binned region (so we don't lose mixed evidence types).
#[derive(Default, Clone, Debug)]
struct EvidenceCounts {
    discordant: usize,
    soft_clip: usize,
    rare_kmer: usize,
}

/// Find candidate breakpoint "hotspots" inside a bin.
/// - discordant pairs → DEL candidates
/// - soft clips → INS / breakpoint candidates
/// - rare kmers → optional fallback (currently disabled by threshold)
pub fn discover_hotspots(reads: &[bam::Record], mindex: &MinimizerIndex, k: usize) -> Vec<Hotspot> {
    let mut evidence_map: HashMap<usize, EvidenceCounts> = HashMap::new();

    // Per-evidence thresholds
    const MIN_RP_SUPPORT: usize = 2;        // discordant pairs: allow low support (important!)
    const MIN_SC_SUPPORT: usize = 1;        // soft-clips: even 1 is informative in synthetic tests
    const MIN_RARE_SUPPORT: usize = 9999;   // keep disabled for now

    // Window sizes by evidence type
    const WIN_DISCORDANT: usize = 500; // key fix: avoids splitting evidence across too many 100bp bins
    const WIN_SOFTCLIP: usize = 100;
    const WIN_RARE: usize = 200;

    println!("--- Starting Hotspot Discovery ---");
    println!("Total reads to process: {}", reads.len());

    for r in reads {
        // keep it simple for now; for real data you may also want r.is_supplementary()
        if r.is_secondary() || r.is_duplicate() {
            continue;
        }

        let pos = r.pos().max(0) as usize;

        // 1) Discordant pairs (DEL signal)
        if r.is_paired() && !r.is_unmapped() && !r.is_mate_unmapped() {
            let isize = r.insert_size().abs() as usize;

            // this condition matches your original intention: abnormal TLEN or weird orientation
            if isize > 1000 || isize < 100 || r.is_reverse() == r.is_mate_reverse() {
                let bin = (pos / WIN_DISCORDANT) * WIN_DISCORDANT;
                evidence_map.entry(bin).or_default().discordant += 1;
                continue; // avoid counting same read as soft-clip/rare-kmer too
            }
        }

        // 2) Soft clips (INS / breakpoint signal)
        {
            let cigar = r.cigar();
            let mut sc_ok = false;

            if let Some(Cigar::SoftClip(n)) = cigar.iter().next() {
                if *n as usize >= k {
                    sc_ok = true;
                }
            }
            if !sc_ok {
                if let Some(Cigar::SoftClip(n)) = cigar.iter().last() {
                    if *n as usize >= k {
                        sc_ok = true;
                    }
                }
            }

            if sc_ok {
                let bin = (pos / WIN_SOFTCLIP) * WIN_SOFTCLIP;
                evidence_map.entry(bin).or_default().soft_clip += 1;
                continue;
            }
        }

        // 3) Rare kmers (disabled by threshold but left here for later)
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
                let bin = (pos / WIN_RARE) * WIN_RARE;
                evidence_map.entry(bin).or_default().rare_kmer += 1;
            }
        }
    }

    println!("Candidate sites found before filtering: {}", evidence_map.len());

    // Emit hotspots per evidence type (one region can produce multiple hotspot types)
    let mut final_hs: Vec<Hotspot> = Vec::new();

    for (bin, c) in evidence_map.into_iter() {
        if c.discordant >= MIN_RP_SUPPORT {
            final_hs.push(Hotspot { local_pos: bin, reason: "discordant".into() });
        }
        if c.soft_clip >= MIN_SC_SUPPORT {
            final_hs.push(Hotspot { local_pos: bin, reason: "soft_clip".into() });
        }
        if c.rare_kmer >= MIN_RARE_SUPPORT {
            final_hs.push(Hotspot { local_pos: bin, reason: "rare_kmer".into() });
        }
    }

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
