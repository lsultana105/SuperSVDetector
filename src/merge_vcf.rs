use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use anyhow::{Context, Result};
use rust_htslib::bcf::{self, header::Header, Writer};
use serde_json::Deserializer;

use crate::confirm::SvCall;
fn add_contigs_from_fai(hdr: &mut Header, ref_fa: &str) -> Result<()> {
    let fai_path = format!("{}.fai", ref_fa);
    let file = File::open(&fai_path)
        .with_context(|| format!("Failed to open FASTA index at {}", fai_path))?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let mut parts = line.split('\t');
        let name = parts.next().unwrap_or_default();
        let len_str = parts.next().unwrap_or_default();
        if name.is_empty() || len_str.is_empty() { continue; }
        let len: u64 = len_str.parse().unwrap_or(0);
        if len > 0 {
            let rec = format!("##contig=<ID={},length={}>", name, len);
            hdr.push_record(rec.as_bytes());
        }
    }
    Ok(())
}

fn main() -> Result<()> {
    let outdir = "demo_results";
    let vcf_path = "demo_results/demo_chr22_standalone.vcf";
    let ref_fa = "dataToTest/demo_ref_chr22.fa";

    // 1. Collect all calls
    let mut calls: Vec<SvCall> = Vec::new();
    for entry in fs::read_dir(outdir)? {
        let path = entry?.path();
        if path.extension().map(|x| x == "json").unwrap_or(false) {
            let file = File::open(&path)?;
            let reader = BufReader::new(file);
            let stream = Deserializer::from_reader(reader).into_iter::<Vec<SvCall>>();
            for v in stream {
                calls.extend(v?);
            }
        }
    }

    // 2. Create header
    let mut hdr = Header::new();
    hdr.push_record(b"##fileformat=VCFv4.2");
    hdr.push_record(b"##source=standalone_merge");
    hdr.push_record(format!("##reference={}", ref_fa).as_bytes());
    hdr.push_record(b"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">");
    hdr.push_record(b"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">");
    hdr.push_record(b"##INFO=<ID=SCORE,Number=1,Type=Integer,Description=\"Score\">");
    hdr.push_record(b"##INFO=<ID=REASON,Number=1,Type=String,Description=\"Hotspot reason\">");

    add_contigs_from_fai(&mut hdr, ref_fa)?;

    hdr.push_record(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");

    // 3. Writer (plain-text VCF)
    let mut w = Writer::from_path(vcf_path, &hdr, false, bcf::Format::Vcf)?;

    // 4. Write records
    for (i, c) in calls.iter().enumerate() {
        let rid = match w.header().name2rid(c.chrom.as_bytes()) {
            Ok(r) => r,
            Err(_) => {
                eprintln!("Warning: chrom '{}' not in header, skipping", c.chrom);
                continue;
            }
        };
        let mut rec = w.empty_record();
        rec.set_rid(Some(rid));
        rec.set_pos((c.pos - 1) as i64);
        rec.set_id(format!("svx_{}", i).as_bytes())?;
        rec.set_qual(0.0);

        let alt = format!("<{}>", c.svtype);
        let alleles = [b"N".as_ref(), alt.as_bytes()];
        rec.set_alleles(&alleles)?;

        rec.push_info_string(b"SVTYPE", &[c.svtype.as_bytes()])?;
        rec.push_info_integer(b"SVLEN", &[c.len as i32])?;
        rec.push_info_integer(b"SCORE", &[c.score])?;
        rec.push_info_string(b"REASON", &[c.reason.as_bytes()])?;
        w.write(&rec)?;
    }

    println!("Standalone VCF written to {}", vcf_path);
    Ok(())
}
