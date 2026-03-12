use anyhow::{anyhow, Context, Result};
use bio::io::fasta;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

use crate::bins::Bin;

pub struct WorkerIo {
    pub bam: bam::IndexedReader,
    pub fasta: fasta::IndexedReader<File>,
    pub tid_cache: HashMap<String, u32>,
    pub tid_names: Vec<String>,
}

impl WorkerIo {
    pub fn new<P: AsRef<Path> + std::fmt::Debug>(bam_path: P, fasta_path: P) -> Result<Self> {
        let bam = bam::IndexedReader::from_path(bam_path.as_ref())
            .with_context(|| format!("failed to open BAM: {}", bam_path.as_ref().display()))?;

        let fasta = fasta::IndexedReader::from_file(&fasta_path)
            .map_err(|e| anyhow!("failed to open FASTA {}: {}", fasta_path.as_ref().display(), e))?;

        let header = bam.header().to_owned();
        let mut tid_names = Vec::new();
        for tid in 0..header.target_count() {
            let name = String::from_utf8_lossy(header.tid2name(tid)).to_string();
            tid_names.push(name);
        }

        Ok(Self {
            bam,
            fasta,
            tid_cache: HashMap::new(),
            tid_names,
        })
    }

    pub fn tid_for_chrom(&mut self, chrom: &str) -> Result<u32> {
        if let Some(&tid) = self.tid_cache.get(chrom) {
            return Ok(tid);
        }

        let tid = self
            .bam
            .header()
            .tid(chrom.as_bytes())
            .ok_or_else(|| anyhow!("chromosome '{}' not found in BAM header", chrom))?;

        self.tid_cache.insert(chrom.to_string(), tid);
        Ok(tid)
    }
}

pub fn fetch_reference_window_with_reader(
    fasta_reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    start: u64,
    end: u64,
) -> Result<Vec<u8>> {
    fasta_reader
        .fetch(chrom, start, end)
        .map_err(|e| anyhow!("FASTA fetch failed for {}:{}-{}: {}", chrom, start, end, e))?;

    let mut seq = Vec::new();
    fasta_reader
        .read(&mut seq)
        .map_err(|e| anyhow!("FASTA read failed for {}:{}-{}: {}", chrom, start, end, e))?;

    Ok(seq)
}

pub fn fetch_reads_overlapping_with_reader(
    io: &mut WorkerIo,
    chrom: &str,
    start: u64,
    end: u64,
) -> Result<Vec<bam::Record>> {
    let tid = io.tid_for_chrom(chrom)?;

    io.bam
        .fetch((tid, start, end))
        .map_err(|e| anyhow!("BAM fetch failed for {}:{}-{}: {}", chrom, start, end, e))?;

    let mut reads = Vec::new();
    for rec_result in io.bam.records() {
        let rec = rec_result.map_err(|e| anyhow!("error reading BAM record: {}", e))?;
        reads.push(rec);
    }

    Ok(reads)
}

pub fn read_bins<P: AsRef<Path>>(bins_path: P) -> Result<Vec<Bin>> {
    let text = std::fs::read_to_string(&bins_path)
        .with_context(|| format!("failed to read bins file: {}", bins_path.as_ref().display()))?;

    let mut bins = Vec::new();

    for (i, line) in text.lines().enumerate() {
        let line_no = i + 1;
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() != 3 {
            return Err(anyhow!(
                "invalid bins line {} in {}: expected 3 columns, got {} -> '{}'",
                line_no,
                bins_path.as_ref().display(),
                parts.len(),
                line
            ));
        }

        let chrom = parts[0].to_string();
        let start: u64 = parts[1]
            .parse()
            .with_context(|| format!("invalid start at line {}: '{}'", line_no, parts[1]))?;
        let end: u64 = parts[2]
            .parse()
            .with_context(|| format!("invalid end at line {}: '{}'", line_no, parts[2]))?;

        if end <= start {
            return Err(anyhow!(
                "invalid bin at line {}: end ({}) must be > start ({})",
                line_no,
                end,
                start
            ));
        }

        bins.push(Bin { chrom, start, end });
    }

    Ok(bins)
}

pub fn estimate_insert_stats<P: AsRef<Path>>(bam_path: P) -> Result<(f64, f64)> {
    let mut reader = bam::Reader::from_path(bam_path.as_ref())
        .with_context(|| format!("failed to open BAM for insert stats: {}", bam_path.as_ref().display()))?;

    const MAX_READS: usize = 200_000;
    const MIN_MAPQ: u8 = 20;

    let mut vals: Vec<f64> = Vec::new();

    for rec_result in reader.records().take(MAX_READS) {
        let rec = rec_result?;
        if rec.is_secondary() || rec.is_supplementary() || rec.is_duplicate() {
            continue;
        }
        if rec.mapq() < MIN_MAPQ {
            continue;
        }
        if !(rec.is_paired() && !rec.is_unmapped() && !rec.is_mate_unmapped()) {
            continue;
        }
        if rec.tid() != rec.mtid() {
            continue;
        }

        let tlen = rec.insert_size().abs();
        if tlen > 0 {
            vals.push(tlen as f64);
        }
    }

    if vals.len() < 10 {
        return Ok((0.0, 0.0));
    }

    vals.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let n = vals.len();
    let lo = n / 20;
    let hi = n - (n / 20);
    let trimmed = &vals[lo..hi];

    if trimmed.is_empty() {
        return Ok((0.0, 0.0));
    }

    let mean = trimmed.iter().sum::<f64>() / trimmed.len() as f64;
    let var = trimmed
        .iter()
        .map(|x| {
            let d = *x - mean;
            d * d
        })
        .sum::<f64>()
        / trimmed.len() as f64;
    let sd = var.sqrt();

    Ok((mean, sd))
}

pub fn get_contig_len_from_fai<P: AsRef<Path>>(ref_fa: P, chrom: &str) -> Result<u64> {
    use std::io::{BufRead, BufReader};

    let fai_path = format!("{}.fai", ref_fa.as_ref().display());
    let file = File::open(&fai_path)
        .with_context(|| format!("failed to open FASTA index {}", fai_path))?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 && parts[0] == chrom {
            return Ok(parts[1].parse::<u64>()?);
        }
    }

    Err(anyhow!(
        "chromosome {} not found in FASTA index {}",
        chrom,
        fai_path
    ))
}