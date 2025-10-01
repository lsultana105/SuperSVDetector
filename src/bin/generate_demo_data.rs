use std::fs::File;
use std::io::{Write, BufWriter};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

// ==== Demo-friendly constants ====
const REF_LENGTH: usize = 2_000_000; // 2 Mb chunk
const LINE_WIDTH: usize = 80;
const READ_LENGTH: usize = 150;
const COVERAGE: usize = 5; // ~5x coverage
const BIN_SIZE: usize = 500_000; // 0.5 Mb per bin
const PAD: usize = 0;
const SVS_PER_BIN: usize = 2; // fewer SVs per bin for speed

fn main() -> std::io::Result<()> {
    std::fs::create_dir_all("dataToTest")?;

    // ================= Reference =================
    let ref_path = "dataToTest/demo_ref_chr22.fa";
    let ref_file = File::create(ref_path)?;
    let mut ref_writer = BufWriter::new(ref_file);

    writeln!(ref_writer, ">chr22_demo_segment")?;

    let mut rng = StdRng::seed_from_u64(42);
    let nucleotides = b"ACGT";
    let mut reference: Vec<char> = Vec::with_capacity(REF_LENGTH);

    for _ in 0..REF_LENGTH {
        reference.push(nucleotides[rng.gen_range(0..4)] as char);
    }

    // Write reference in FASTA format
    for (i, &nt) in reference.iter().enumerate() {
        write!(ref_writer, "{}", nt)?;
        if (i + 1) % LINE_WIDTH == 0 {
            writeln!(ref_writer)?;
        }
    }
    writeln!(ref_writer)?; // final newline
    println!("Reference written to {}", ref_path);

    // ================= Bins =================
    let bins_path = "dataToTest/demo_bins_chr22.txt";
    let mut bins_file = BufWriter::new(File::create(bins_path)?);

    let mut start = 0;
    while start < REF_LENGTH {
        let end = std::cmp::min(start + BIN_SIZE, REF_LENGTH);
        writeln!(bins_file, "chr22_demo_segment\t{}\t{}", start, end)?;
        start = end;
    }
    println!("Bins written to {}", bins_path);

    // ================= SV positions =================
    let mut sv_positions = Vec::new();
    let num_bins = ((REF_LENGTH + BIN_SIZE - 1) / BIN_SIZE) as usize;

    for bin_idx in 0..num_bins {
        let bin_start = bin_idx * BIN_SIZE;
        for _ in 0..SVS_PER_BIN {
            let pos = bin_start + rng.gen_range(0..BIN_SIZE.min(REF_LENGTH - bin_start));
            let len = rng.gen_range(20..80); // shorter SVs for speed
            sv_positions.push((pos, len));
        }
    }
    println!("Generated {} SV positions across bins", sv_positions.len());

    // ================= Reads with simulated SVs =================
    let reads_path = "dataToTest/demo_reads_chr22.fastq";
    let reads_file = File::create(reads_path)?;
    let mut reads_writer = BufWriter::new(reads_file);

    let num_reads = REF_LENGTH * COVERAGE / READ_LENGTH;
    for i in 0..num_reads {
        let start = rng.gen_range(0..(REF_LENGTH - READ_LENGTH));
        let mut seq: Vec<char> = reference[start..start + READ_LENGTH].to_vec();

        // With small probability, introduce a structural variant
        let roll: f64 = rng.gen();
        if roll < 0.02 {
            // deletion
            let del_len = rng.gen_range(10..30).min(seq.len() / 2);
            let del_pos = rng.gen_range(0..seq.len() - del_len);
            seq.splice(del_pos..del_pos + del_len, std::iter::empty());
        } else if roll < 0.04 {
            // insertion
            let ins_len = rng.gen_range(10..30);
            let ins_pos = rng.gen_range(0..seq.len());
            let ins_seq: String = (0..ins_len)
                .map(|_| nucleotides[rng.gen_range(0..4)] as char)
                .collect();
            seq.splice(ins_pos..ins_pos, ins_seq.chars());
        } else if roll < 0.06 {
            // inversion
            let inv_len = rng.gen_range(10..30).min(seq.len() / 2);
            let inv_pos = rng.gen_range(0..seq.len() - inv_len);
            let mut segment: Vec<char> = seq[inv_pos..inv_pos + inv_len].to_vec();
            segment.reverse();
            seq.splice(inv_pos..inv_pos + inv_len, segment);
        }

        let seq_str: String = seq.iter().collect();

        writeln!(reads_writer, "@read{}", i)?;
        writeln!(reads_writer, "{}", seq_str)?;
        writeln!(reads_writer, "+")?;
        writeln!(reads_writer, "{}", "I".repeat(seq_str.len()))?;
    }
    println!("FASTQ reads with simulated SVs written to {}", reads_path);

    // ================= Hints for next steps =================
    println!("\nNow run:");
    println!("samtools faidx {}", ref_path);
    println!("bwa index {}", ref_path);
    println!("bwa mem {} {} | samtools view -Sb - | samtools sort -o dataToTest/demo_reads_chr22.sorted.bam", ref_path, reads_path);
    println!("samtools index dataToTest/demo_reads_chr22.sorted.bam");
    println!("mpirun -np 2 target/release/SuperSVDetector call --ref-fa {} --bam dataToTest/demo_reads_chr22.sorted.bam --bins {} --outdir demo_results --mpi", ref_path, bins_path);

    Ok(())
}
