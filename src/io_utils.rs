use anyhow::{Result, anyhow};
use rust_htslib::bam::{self, Read};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use bio::io::fasta;

use crate::bins::Bin;

pub fn tid_to_name(reader: &bam::IndexedReader, tid: i32) -> Option<String> {
    if tid < 0 { return None; }
    let hdr = reader.header().to_owned();
    let bytes: &[u8] = hdr.tid2name(tid as u32);
    Some(String::from_utf8_lossy(bytes).to_string())
}

/// Fetch a reference subsequence [start, end) for a given chromosome.
pub fn fetch_reference_window<P: AsRef<Path> + std::fmt::Debug>(
    ref_fa: P,
    chrom: &str,
    start: u64,
    end: u64,
) -> Result<Vec<u8>> {
    let mut rdr = fasta::IndexedReader::from_file(&ref_fa)?;

    rdr.fetch(chrom, start, end)?;

    let mut seq: Vec<u8> = Vec::new();
    rdr.read(&mut seq)?;

    for b in seq.iter_mut() {
        *b = b.to_ascii_uppercase();
    }

    Ok(seq)
}

pub fn open_bam(path: &str) -> Result<bam::IndexedReader> {
    Ok(bam::IndexedReader::from_path(path)?)
}

pub fn fetch_reads_overlapping(
    reader: &mut bam::IndexedReader,
    chrom: &str,
    start: u64,
    end: u64,
) -> Result<Vec<bam::Record>> {
    let header = reader.header().to_owned();
    let tid = header
        .tid(chrom.as_bytes())
        .ok_or_else(|| anyhow!("chrom not found: {chrom}"))? as u32;

    reader.fetch((tid, start as i64, end as i64))?;
    let mut v = Vec::new();
    for r in reader.records() {
        v.push(r?);
    }
    Ok(v)
}

pub fn read_bins(bins_path: &str) -> Result<Vec<Bin>> {
    let fh = File::open(bins_path)?;
    let br = BufReader::new(fh);
    let mut out = Vec::new();
    for line in br.lines() {
        let l = line?;
        if l.trim().is_empty() {
            continue;
        }
        let parts: Vec<_> = l.split('\t').collect();
        let b = Bin {
            chrom: parts[0].to_string(),
            start: parts[1].parse()?,
            end: parts[2].parse()?,
        };
        out.push(b);
    }
    Ok(out)
}


/// Estimate insert size mean and standard deviation (TLEN) for proper pairs.
/// Filters extreme TLENs to avoid insane values.
///
/// Returns (mean, sd). If insufficient data, returns (0.0, 0.0).
pub fn estimate_insert_stats(bam_path: &str) -> Result<(f64, f64)> {
    let mut reader = bam::Reader::from_path(bam_path)?;

    let mut vals: Vec<f64> = Vec::new();
    vals.reserve(200_000);

    for rec in reader.records() {
        let r = rec?;

        // keep only clean-ish proper pairs
        if !r.is_paired() || r.is_unmapped() || r.is_mate_unmapped() {
            continue;
        }
        if r.is_secondary() || r.is_duplicate() {
            continue;
        }
        if r.tid() != r.mtid() {
            continue;
        }
        if !r.is_proper_pair() {
            continue;
        }

        let tlen = r.insert_size().abs() as f64;

        // Keep plausible values (tune if needed)
        if tlen <= 0.0 || tlen > 2000.0 {
            continue;
        }

        vals.push(tlen);

        // cap sample size for speed (optional)
        if vals.len() >= 200_000 {
            break;
        }
    }

    if vals.len() < 2000 {
        return Ok((0.0, 0.0));
    }

    let n = vals.len() as f64;
    let mean = vals.iter().sum::<f64>() / n;

    let var = vals
        .iter()
        .map(|x| {
            let d = x - mean;
            d * d
        })
        .sum::<f64>()
        / (n - 1.0);

    let sd = var.sqrt();

    Ok((mean, sd))
}

