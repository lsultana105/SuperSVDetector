use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use anyhow::{bail, Context, Result};
use rust_htslib::bcf::{self, header::Header, Writer};
use serde_json::Deserializer;
use crate::cluster::cluster_calls;
use crate::confirm::SvCall;

pub fn write_bin_json(outdir: &str, bin: &crate::bins::Bin, calls: &[SvCall]) -> Result<()> {
    fs::create_dir_all(outdir)?;
    let path = format!("{outdir}/{}_{}_{}.json", bin.chrom, bin.start, bin.end);
    let mut f = BufWriter::new(File::create(path)?);
    let s = serde_json::to_string(&calls)?;
    f.write_all(s.as_bytes())?;
    Ok(())
}

fn add_contigs_from_fai(hdr: &mut Header, ref_fa: &str) -> Result<()> {
    let fai_path = format!("{}.fai", ref_fa);
    let file = File::open(&fai_path).with_context(|| {
        format!(
            "Failed to open FASTA index ({}.fai). Did you run `samtools faidx {ref_fa}`?",
            ref_fa
        )
    })?;
    let reader = BufReader::new(file);
    let mut n_contigs = 0usize;

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let mut it = line.split('\t');
        let name = it.next().unwrap_or_default();
        let len_str = it.next().unwrap_or_default();
        if name.is_empty() || len_str.is_empty() {
            continue;
        }

        let len: u64 = len_str.parse().unwrap_or(0);
        if len == 0 {
            continue;
        }

        let rec = format!("##contig=<ID={},length={}>", name, len);
        hdr.push_record(rec.as_bytes());
        n_contigs += 1;
    }

    if n_contigs == 0 {
        bail!("No contigs found in {}.fai; header would be missing chromosomes.", ref_fa);
    }
    Ok(())
}

pub fn merge_json_to_vcf(outdir: &str, vcf_path: &str, ref_fa: &str) -> Result<()> {
    // 1) Collect calls
    let mut calls: Vec<SvCall> = Vec::new();
    let mut n_json = 0usize;

    for entry in fs::read_dir(outdir)? {
        let path = entry?.path();
        if path.extension().map(|x| x == "json").unwrap_or(false) {
            let file = File::open(&path)
                .with_context(|| format!("Opening per-bin JSON {}", path.display()))?;
            let reader = BufReader::new(file);
            let stream = Deserializer::from_reader(reader).into_iter::<Vec<SvCall>>();
            for v in stream {
                let v = v.with_context(|| format!("Decoding JSON file {}", path.display()))?;
                calls.extend(v);
            }
            n_json += 1;
        }
    }

    // Cluster nearby calls
    let calls = cluster_calls(calls, 500);

    // 2) Header
    let mut hdr = Header::new();
    hdr.push_record(b"##fileformat=VCFv4.2");
    hdr.push_record(b"##source=SuperSVDetector");
    hdr.push_record(format!("##reference={}", ref_fa).as_bytes());
    hdr.push_record(b"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">");
    hdr.push_record(b"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV (negative for deletions)\">");
    hdr.push_record(b"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of the variant\">");
    hdr.push_record(b"##INFO=<ID=SCORE,Number=1,Type=Integer,Description=\"Support/heuristic score\">");
    hdr.push_record(b"##INFO=<ID=REASON,Number=1,Type=String,Description=\"Evidence source\">");
    // for BND debug/analysis
    hdr.push_record(b"##INFO=<ID=MATECHR,Number=1,Type=String,Description=\"Mate chromosome for BND\">");
    hdr.push_record(b"##INFO=<ID=MATEPOS,Number=1,Type=Integer,Description=\"Mate position (1-based) for BND\">");
    add_contigs_from_fai(&mut hdr, ref_fa)?;
    hdr.push_record(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");

    let mut w = Writer::from_path(vcf_path, &hdr, true, bcf::Format::Vcf)
        .with_context(|| format!("Opening VCF for write at {}", vcf_path))?;

    // 3) Write calls
    let mut written = 0usize;
    for (i, c) in calls.iter().enumerate() {
        let rid = match w.header().name2rid(c.chrom.as_bytes()) {
            Ok(r) => r,
            Err(_) => {
                eprintln!(
                    "Warning: chrom '{}' not found in VCF header; skipping record at {}",
                    c.chrom, c.pos
                );
                continue;
            }
        };

        let mut rec = w.empty_record();
        rec.set_rid(Some(rid));
        rec.set_pos((c.pos.saturating_sub(1)) as i64);
        rec.set_id(format!("svx_{}", i).as_bytes())?;
        rec.set_qual(0.0);

        // ALT formatting
        // - DEL/INS: symbolic
        // - BND: minimal breakend string
        let alt_string = if c.svtype == "BND" {
            if let (Some(mchr), Some(mpos)) = (&c.mate_chrom, c.mate_pos) {
                // minimal: N[chr:pos[
                format!("N[{}:{}[", mchr, mpos)
            } else {
                "<BND>".to_string()
            }
        } else {
            format!("<{}>", c.svtype)
        };

        let alleles = [b"N".as_ref(), alt_string.as_bytes()];
        rec.set_alleles(&alleles)?;

        // INFO
        rec.push_info_string(b"SVTYPE", &[c.svtype.as_bytes()])?;
        rec.push_info_integer(b"SVLEN", &[c.len as i32])?;

        if c.end > i32::MAX as u64 {
            eprintln!("Warning: END too large for i32 ({}), skipping", c.end);
            continue;
        }
        rec.push_info_integer(b"END", &[c.end as i32])?;
        rec.push_info_integer(b"SCORE", &[c.score])?;
        rec.push_info_string(b"REASON", &[c.reason.as_bytes()])?;

        if let Some(mchr) = &c.mate_chrom {
            rec.push_info_string(b"MATECHR", &[mchr.as_bytes()])?;
        }
        if let Some(mpos) = c.mate_pos {
            if mpos <= i32::MAX as u64 {
                rec.push_info_integer(b"MATEPOS", &[mpos as i32])?;
            }
        }

        w.write(&rec)?;
        written += 1;
    }

    println!(
        "Merged {} JSON file(s), wrote {} record(s) to {}",
        n_json, written, vcf_path
    );

    Ok(())
}
