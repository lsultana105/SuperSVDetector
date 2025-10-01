//! Lightweight minimizer index and hotspot discovery.
//! The goal is to quickly highlight candidate regions, not to be perfect.

use hashbrown::HashMap;
use rust_htslib::bam;
use serde::{Serialize, Deserialize};

/// Minimizer index: maps minimizer (u64) → list of positions in reference window
pub struct MinimizerIndex {
    k: usize,
    w: usize,
    pub map: HashMap<u64, Vec<usize>>, // positions within ref_window
}

impl MinimizerIndex {
    /// Build a simple minimizer table over a reference window.
    pub fn build(ref_window: &[u8], k: usize, w: usize) -> Self {
        assert!(w >= k);
        let mut map: HashMap<u64, Vec<usize>> = HashMap::new();
        if ref_window.len() < w {
            return Self { k, w, map };
        }

        // encode A/C/G/T to 2 bits; ignore Ns by skipping
        let enc = |b: u8| match b {
            b'A' => Some(0u64),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        };

        // rolling hash for k-mers inside each window of size w
        for i in 0..=ref_window.len() - w {
            let mut best: Option<(u64, usize)> = None; // (value, start_pos)
            for j in i..=i + w - k {
                // compute k-mer
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
                if !ok {
                    continue;
                }
                if best.map(|(v, _)| val < v).unwrap_or(true) {
                    best = Some((val, j));
                }
            }
            if let Some((val, jpos)) = best {
                map.entry(val).or_default().push(jpos);
            }
        }
        Self { k, w, map }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Hotspot {
    /// genomic coordinate *within the bin window* (offset from bin.start)
    pub local_pos: usize,
    /// reason for discovery (softclip/minimizer rarity/discordant pairs/etc.)
    pub reason: String,
}

/// Heuristic discovery of candidate SV breakpoints
pub fn discover_hotspots(
    reads: &[bam::Record],
    mindex: &MinimizerIndex,
    _ref_window: &[u8],
    k: usize,
) -> Vec<Hotspot> {
    use rust_htslib::bam::record::Cigar;
    let mut hs = Vec::new();

    // ---- 1) Soft-clips
    for r in reads {
        let mut left_soft = 0usize;
        let mut right_soft = 0usize;
        for (idx, c) in r.cigar().iter().enumerate() {
            if let Cigar::SoftClip(n) = c {
                if idx == 0 {
                    left_soft = *n as usize;
                }
            }
        }
        if let Some(Cigar::SoftClip(n)) = r.cigar().last() {
            right_soft = *n as usize;
        }

        let pos = r.pos().max(0) as usize;
        if left_soft >= k {
            hs.push(Hotspot {
                local_pos: pos,
                reason: format!("left_softclip>={k}"),
            });
        }
        if right_soft >= k {
            hs.push(Hotspot {
                local_pos: pos.saturating_add(r.seq_len() as usize),
                reason: format!("right_softclip>={k}"),
            });
        }
    }

    // ---- 2) Minimizer rarity
    for r in reads {
        let seq = r.seq();
        if seq.len() < k {
            continue;
        }
        let mut rare_count = 0usize;
        let mut total = 0usize;
        for i in 0..=seq.len() - k {
            let val = encode_kmer(&seq, i, k);
            if let Some(v) = val {
                total += 1;
                if mindex.map.get(&v).map(|p| p.len()).unwrap_or(0) == 0 {
                    rare_count += 1;
                }
            }
        }
        if total >= 10 && rare_count * 2 > total {
            // >50% rare
            let pos = r.pos().max(0) as usize;
            hs.push(Hotspot {
                local_pos: pos,
                reason: "minimizer_rare".into(),
            });
        }
    }

    // ---- 3) Discordant pairs
    for r in reads {
        if r.is_paired() && !r.is_unmapped() && !r.is_mate_unmapped() {
            let insert_size = (r.insert_size().abs()) as usize;

            // Adjust thresholds as needed
            if insert_size > 1000 {
                hs.push(Hotspot {
                    local_pos: r.pos().max(0) as usize,
                    reason: format!("discordant_large_insert={}", insert_size),
                });
            }
            if insert_size < 100 {
                hs.push(Hotspot {
                    local_pos: r.pos().max(0) as usize,
                    reason: format!("discordant_small_insert={}", insert_size),
                });
            }

            // Orientation check: normally FR; abnormal = candidate inversion/translocation
            if r.is_reverse() == r.is_mate_reverse() {
                hs.push(Hotspot {
                    local_pos: r.pos().max(0) as usize,
                    reason: "discordant_orientation".into(),
                });
            }
        }
    }

    // ---- Deduplicate coarse positions
    hs.sort_by_key(|h| h.local_pos);
    hs.dedup_by(|a, b| (a.local_pos / 50) == (b.local_pos / 50));
    hs
}

#[inline]
fn encode_kmer(seq: &rust_htslib::bam::record::Seq, i: usize, k: usize) -> Option<u64> {
    let mut val = 0u64;
    for t in 0..k {
        val = (val << 2)
            | match nuc_to_acgt(seq.as_bytes()[(i + t) as usize]) {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
    }
    Some(val)
}

#[inline]
fn nuc_to_acgt(n: u8) -> u8 {
    match n {
        1 => b'A',
        2 => b'C',
        4 => b'G',
        8 => b'T',
        _ => b'N',
    }
}
