use anyhow::Result;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Bin {
    pub chrom: String,
    pub start: u64,
    pub end: u64
}

pub fn build_bins<P: AsRef<Path> + std::fmt::Debug>(ref_fa: P, bin: usize, _pad: usize, out: &str) -> Result<()> {
    // Load the FASTA index (.fai)
    let fai_path = format!("{}.fai", ref_fa.as_ref().display());
    let fai = bio::io::fasta::Index::from_file(&fai_path)?;
    let mut writer = BufWriter::new(File::create(out)?);

    for record in fai.sequences() {
        let name = record.name.clone();
        let len = record.len as usize;
        let mut s = 0usize;
        while s < len {
            let e = (s + bin).min(len);
            writeln!(writer, "{name}\t{s}\t{e}")?;
            s = e;
        }
    }
    Ok(())
}
