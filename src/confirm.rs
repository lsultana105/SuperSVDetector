use anyhow::Result;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

use crate::hotspot::Hotspot;
use crate::suffix::SuffixWorkspace;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SvCall {
    pub chrom: String,
    pub pos: u64,
    pub end: u64,
    pub svtype: String,
    pub len: isize,
    pub score: i32,
    pub reason: String,

    // Minimal BND support (optional)
    pub mate_chrom: Option<String>,
    pub mate_pos: Option<u64>,
}

// median for u64
fn median_u64(v: &mut Vec<u64>) -> u64 {
    v.sort_unstable();
    v[v.len() / 2]
}
fn median_isize(v: &mut Vec<isize>) -> isize {
    v.sort_unstable();
    v[v.len() / 2]
}

pub fn confirm_breakpoint(
    bin_start: u64,
    chrom: &str,
    _ref_window: &[u8],
    _sfx: &SuffixWorkspace,
    hs: &Hotspot,
    band: usize,
    insert_mean: f64,
    insert_sd: f64,
    reads: &[bam::Record],
    k: usize,
    tid_names: &[String], // tid -> chrom name
) -> Result<Option<SvCall>> {
    // Hotspot local -> global window
    let center = hs.local_pos;
    let local_start = center.saturating_sub(band);
    let local_end = center + band;

    let g_start0 = bin_start + local_start as u64;
    let g_end0 = bin_start + local_end as u64;

    // Gather overlapping reads
    let mut support: Vec<&bam::Record> = Vec::new();
    for r in reads.iter() {
        let pos0 = r.pos().max(0) as u64;
        let end0 = pos0 + r.seq_len() as u64; // coarse gate
        if end0 > g_start0 && pos0 < g_end0 {
            support.push(r);
        }
    }
    if support.is_empty() {
        return Ok(None);
    }

    const MIN_MAPQ: u8 = 20;

    // -------------------------
    // helper: reference consumed by cigar (M,=,X,D,N)
    // -------------------------
    fn ref_consumed(cig: &bam::record::CigarStringView) -> u64 {
        use bam::record::Cigar::*;
        let mut n: u64 = 0;
        for c in cig.iter() {
            match *c {
                Match(l) | Equal(l) | Diff(l) | Del(l) | RefSkip(l) => n += l as u64,
                Ins(_) | SoftClip(_) | HardClip(_) | Pad(_) => {}
            }
        }
        n
    }

    // -------------------------
    // helper: parse first SA entry "chr,pos,strand,cigar,mapq,nm;..."
    // -------------------------
    fn parse_sa_first(sa: &str) -> Option<(String, u64, String)> {
        let first = sa.split(';').next()?.trim();
        if first.is_empty() {
            return None;
        }
        let mut it = first.split(',');
        let chr = it.next()?.to_string();
        let pos1: u64 = it.next()?.parse().ok()?; // 1-based
        let _strand = it.next()?; // not used
        let cigar = it.next()?.to_string();
        Some((chr, pos1, cigar))
    }

    // -------------------------
    // helper: parse SA cigar ref-consumed (quick parser; avoids CigarString::try_from)
    // We only need how much reference it consumes.
    // -------------------------
    fn sa_ref_consumed(cigar: &str) -> Option<u64> {
        let mut n: u64 = 0;
        let mut num: u64 = 0;
        let mut any = false;

        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                num = num * 10 + (ch as u8 - b'0') as u64;
                any = true;
            } else {
                if !any {
                    return None;
                }
                match ch {
                    'M' | '=' | 'X' | 'D' | 'N' => n += num,
                    'I' | 'S' | 'H' | 'P' => {}
                    _ => return None,
                }
                num = 0;
                any = false;
            }
        }
        // cigar should end in an op; if not, fail
        if any { return None; }
        Some(n)
    }

    // -------------------------
    // DEL from discordant PE
    // -------------------------
    if hs.reason == "discordant_del" {
        const MIN_DEL_BP: u64 = 50;
        const MIN_PE_SUPPORT: usize = 6;
        const Z_CUTOFF: f64 = 5.0;
        const MAX_BP_SPREAD: u64 = 300;

        if insert_sd <= 0.0 {
            return Ok(None);
        }

        let mut left_bps: Vec<u64> = Vec::new();
        let mut right_bps: Vec<u64> = Vec::new();
        let mut seen_qname: HashSet<Vec<u8>> = HashSet::new();

        for r in &support {
            if r.is_secondary() || r.is_duplicate() {
                continue;
            }
            if r.mapq() < MIN_MAPQ {
                continue;
            }
            if !(r.is_paired() && !r.is_unmapped() && !r.is_mate_unmapped()) {
                continue;
            }
            if r.tid() != r.mtid() {
                continue;
            }
            if r.is_proper_pair() {
                continue;
            }

            // one per pair
            if r.pos() >= r.mpos() {
                continue;
            }

            let tlen = r.insert_size().abs() as f64;
            let z = (tlen - insert_mean) / insert_sd;
            if z < Z_CUTOFF {
                continue;
            }

            let qn = r.qname().to_vec();
            if !seen_qname.insert(qn) {
                continue;
            }

            let rpos0 = r.pos().max(0) as u64;
            let mpos0 = r.mpos().max(0) as u64;

            let left0 = rpos0.min(mpos0);
            let right0 = rpos0.max(mpos0);

            let read_len = r.seq_len() as u64;
            let left_bp0 = left0 + read_len;
            let right_bp0 = right0;

            if right_bp0 > left_bp0 {
                left_bps.push(left_bp0);
                right_bps.push(right_bp0);
            }
        }

        if left_bps.len() < MIN_PE_SUPPORT || right_bps.len() < MIN_PE_SUPPORT {
            return Ok(None);
        }

        left_bps.sort_unstable();
        right_bps.sort_unstable();
        let left_spread = left_bps[left_bps.len() - 1] - left_bps[0];
        let right_spread = right_bps[right_bps.len() - 1] - right_bps[0];
        if left_spread > MAX_BP_SPREAD || right_spread > MAX_BP_SPREAD {
            return Ok(None);
        }

        let start0 = median_u64(&mut left_bps);
        let end0 = median_u64(&mut right_bps);

        if end0 <= start0 || (end0 - start0) < MIN_DEL_BP {
            return Ok(None);
        }

        let svlen = -((end0 - start0) as isize);
        return Ok(Some(SvCall {
            chrom: chrom.into(),
            pos: start0 + 1,
            end: end0 + 1,
            svtype: "DEL".into(),
            len: svlen,
            score: left_bps.len() as i32,
            reason: "PE".into(),
            mate_chrom: None,
            mate_pos: None,
        }));
    }

    // -------------------------
    // BND from inter-chrom PE (minimal)
    // -------------------------
    if hs.reason == "bnd_pe" {
        const MIN_BND_SUPPORT: usize = 3;
        const MAX_POS_SPREAD: u64 = 500;

        let mut a_pos: Vec<u64> = Vec::new();
        let mut b_pos: Vec<u64> = Vec::new();
        let mut mate_chrom_name: Option<String> = None;
        let mut seen_qname: HashSet<Vec<u8>> = HashSet::new();

        for r in &support {
            if r.is_secondary() || r.is_duplicate() {
                continue;
            }
            if r.mapq() < MIN_MAPQ {
                continue;
            }
            if !(r.is_paired() && !r.is_unmapped() && !r.is_mate_unmapped()) {
                continue;
            }
            if r.tid() == r.mtid() {
                continue;
            }

            // one per pair
            if r.pos() >= r.mpos() {
                continue;
            }

            let qn = r.qname().to_vec();
            if !seen_qname.insert(qn) {
                continue;
            }

            let rpos0 = r.pos().max(0) as u64;
            let mpos0 = r.mpos().max(0) as u64;
            a_pos.push(rpos0);
            b_pos.push(mpos0);

            let mtid = r.mtid();
            if mtid >= 0 {
                let mtid = mtid as usize;
                if mtid < tid_names.len() {
                    mate_chrom_name = Some(tid_names[mtid].clone());
                }
            }
        }

        if a_pos.len() < MIN_BND_SUPPORT || mate_chrom_name.is_none() {
            return Ok(None);
        }

        a_pos.sort_unstable();
        b_pos.sort_unstable();
        let a_spread = a_pos[a_pos.len() - 1] - a_pos[0];
        let b_spread = b_pos[b_pos.len() - 1] - b_pos[0];
        if a_spread > MAX_POS_SPREAD || b_spread > MAX_POS_SPREAD {
            return Ok(None);
        }

        let pos0 = median_u64(&mut a_pos);
        let mate0 = median_u64(&mut b_pos);

        return Ok(Some(SvCall {
            chrom: chrom.into(),
            pos: pos0 + 1,
            end: pos0 + 1,
            svtype: "BND".into(),
            len: 0,
            score: a_pos.len() as i32,
            reason: "PE".into(),
            mate_chrom: mate_chrom_name,
            mate_pos: Some(mate0 + 1),
        }));
    }

    // -------------------------
    // SR hotspot: try DEL-from-SA first, else INS-from-softclip
    // -------------------------
    if hs.reason == "soft_clip_sr" {
        // DEL-from-SA
        const MIN_DEL_BP: u64 = 50;
        const MIN_SR_DEL_SUPPORT: usize = 4;
        const MAX_BP_SPREAD_DEL: u64 = 300;

        // INS-from-softclip (fallback)
        const MIN_SR_SUPPORT: usize = 4;
        const MAX_POS_SPREAD_INS: u64 = 150;
        const MIN_INS_LEN: isize = 50;

        let min_clip_sr: u32 = (k as u32).max(20);

        // -------------------------
        // (1) DEL from SA split reads
        // -------------------------
        {
            let mut left_bps: Vec<u64> = Vec::new();
            let mut right_bps: Vec<u64> = Vec::new();
            let mut seen_qname: HashSet<Vec<u8>> = HashSet::new();

            for r in &support {
                if r.is_secondary() || r.is_duplicate() {
                    continue;
                }
                if r.mapq() < MIN_MAPQ {
                    continue;
                }

                // ✅ FIX: Aux handling (no .string() in your rust_htslib)
                let sa_str: &str = match r.aux(b"SA") {
                    Ok(bam::record::Aux::String(s)) => s,
                    _ => continue,
                };

                let qn = r.qname().to_vec();
                if !seen_qname.insert(qn) {
                    continue;
                }

                let (sa_chr, sa_pos1, sa_cigar_str) = match parse_sa_first(sa_str) {
                    Some(x) => x,
                    None => continue,
                };
                if sa_chr != chrom {
                    continue;
                }

                // primary alignment ref interval
                let rpos0 = r.pos().max(0) as u64;
                let r_end0 = rpos0 + ref_consumed(&r.cigar());

                // SA alignment ref interval
                let sa_pos0 = sa_pos1.saturating_sub(1);
                let sa_ref = match sa_ref_consumed(&sa_cigar_str) {
                    Some(x) => x,
                    None => continue,
                };
                let sa_end0 = sa_pos0 + sa_ref;

                // deletion gap between the two blocks
                let (a0, a1) = (rpos0, r_end0);
                let (b0, b1) = (sa_pos0, sa_end0);

                let (left_end, right_start) = if a0 <= b0 { (a1, b0) } else { (b1, a0) };

                if right_start <= left_end {
                    continue;
                }
                let del_len = right_start - left_end;
                if del_len < MIN_DEL_BP {
                    continue;
                }

                left_bps.push(left_end);
                right_bps.push(right_start);
            }

            if left_bps.len() >= MIN_SR_DEL_SUPPORT && right_bps.len() >= MIN_SR_DEL_SUPPORT {
                left_bps.sort_unstable();
                right_bps.sort_unstable();

                let l_spread = left_bps[left_bps.len() - 1] - left_bps[0];
                let r_spread = right_bps[right_bps.len() - 1] - right_bps[0];

                if l_spread <= MAX_BP_SPREAD_DEL && r_spread <= MAX_BP_SPREAD_DEL {
                    let start0 = median_u64(&mut left_bps);
                    let end0 = median_u64(&mut right_bps);

                    if end0 > start0 && (end0 - start0) >= MIN_DEL_BP {
                        let svlen = -((end0 - start0) as isize);
                        return Ok(Some(SvCall {
                            chrom: chrom.into(),
                            pos: start0 + 1,
                            end: end0 + 1,
                            svtype: "DEL".into(),
                            len: svlen,
                            score: left_bps.len() as i32,
                            reason: "SR".into(),
                            mate_chrom: None,
                            mate_pos: None,
                        }));
                    }
                }
            }
        }

        // -------------------------
        // (2) INS from SR-like soft clips (fallback)
        // -------------------------
        let mut clip_lens: Vec<isize> = Vec::new();
        let mut clip_pos0: Vec<u64> = Vec::new();
        let mut seen_qname: HashSet<Vec<u8>> = HashSet::new();

        for r in &support {
            if r.is_secondary() || r.is_duplicate() {
                continue;
            }
            if r.mapq() < MIN_MAPQ {
                continue;
            }

            let split_like = r.is_supplementary() || r.aux(b"SA").is_ok();
            if !split_like {
                continue;
            }

            let qn = r.qname().to_vec();
            if !seen_qname.insert(qn) {
                continue;
            }

            let pos0 = r.pos().max(0) as u64;
            let cig = r.cigar();

            let mut ok = false;

            if let Some(Cigar::SoftClip(n)) = cig.iter().next() {
                if *n >= min_clip_sr {
                    ok = true;
                    clip_lens.push(*n as isize);
                    clip_pos0.push(pos0);
                }
            }
            if !ok {
                if let Some(Cigar::SoftClip(n)) = cig.iter().last() {
                    if *n >= min_clip_sr {
                        ok = true;
                        clip_lens.push(*n as isize);
                        clip_pos0.push(pos0 + r.seq_len() as u64);
                    }
                }
            }

            if ok {
                continue;
            }
        }

        if clip_lens.len() < MIN_SR_SUPPORT {
            return Ok(None);
        }

        clip_pos0.sort_unstable();
        let spread = clip_pos0[clip_pos0.len() - 1] - clip_pos0[0];
        if spread > MAX_POS_SPREAD_INS {
            return Ok(None);
        }

        let approx_len = median_isize(&mut clip_lens);
        if approx_len < MIN_INS_LEN {
            return Ok(None);
        }

        let ins0 = median_u64(&mut clip_pos0);
        let pos_vcf = ins0 + 1;

        return Ok(Some(SvCall {
            chrom: chrom.into(),
            pos: pos_vcf,
            end: pos_vcf,
            svtype: "INS".into(),
            len: approx_len,
            score: clip_lens.len() as i32,
            reason: "SR".into(),
            mate_chrom: None,
            mate_pos: None,
        }));
    }

    Ok(None)
}



fn cigar_ref_consumed(cig: &bam::record::CigarStringView) -> u64 {
    use bam::record::Cigar::*;
    let mut n: u64 = 0;
    for c in cig.iter() {
        match *c {
            Match(l) | Equal(l) | Diff(l) | Del(l) | RefSkip(l) => n += l as u64,
            Ins(_) | SoftClip(_) | HardClip(_) | Pad(_) => {}
        }
    }
    n
}

// Parse first SA entry: "chr,pos,strand,cigar,mapq,nm;..."
fn parse_sa_first(sa: &str) -> Option<(String, u64, char, String)> {
    let first = sa.split(';').next()?.trim();
    if first.is_empty() { return None; }
    let mut it = first.split(',');
    let chr = it.next()?.to_string();
    let pos1: u64 = it.next()?.parse().ok()?; // 1-based
    let strand: char = it.next()?.chars().next()?;
    let cigar = it.next()?.to_string();
    Some((chr, pos1, strand, cigar))
}
