use std::fs::File;
use std::io::{Write, BufWriter};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

// ==== Constants for a 10-20GB BAM Test ====
const REF_LENGTH: usize = 250_000_000; // 250 Mb (Large Human Chromosome size)
const LINE_WIDTH: usize = 80;
const READ_LENGTH: usize = 150;
const COVERAGE: usize = 30;           // High depth for realistic testing
const BIN_SIZE: usize = 1_000_000;     // 1 Mb bins
const INSERT_SIZE: usize = 300;        // Distance between paired reads
const INSERT_STD: usize = 50;         // Variation in distance

fn main() -> std::io::Result<()> {
    let out_dir = "dataToTest";
    std::fs::create_dir_all(out_dir)?;

    // ================= 1. Generate Reference =================
    let ref_path = format!("{}/demo_ref.fa", out_dir);
    let mut ref_writer = BufWriter::new(File::create(&ref_path)?);
    writeln!(ref_writer, ">chr_sim_1")?;

    let mut rng = StdRng::seed_from_u64(42);
    let nucleotides = b"ACGT";
    let mut reference: Vec<char> = Vec::with_capacity(REF_LENGTH);

    for _ in 0..REF_LENGTH {
        reference.push(nucleotides[rng.gen_range(0..4)] as char);
    }

    for (i, &nt) in reference.iter().enumerate() {
        write!(ref_writer, "{}", nt)?;
        if (i + 1) % LINE_WIDTH == 0 { writeln!(ref_writer)?; }
    }
    writeln!(ref_writer)?;
    println!("Reference (250MB) written to {}", ref_path);

    // ================= 2. Generate Bins =================
    let bins_path = format!("{}/demo_bins.txt", out_dir);
    let mut bins_file = BufWriter::new(File::create(bins_path)?);
    let mut start = 0;
    while start < REF_LENGTH {
        let end = std::cmp::min(start + BIN_SIZE, REF_LENGTH);
        writeln!(bins_file, "chr_sim_1\t{}\t{}", start, end)?;
        start = end;
    }

    // ================= 3. Simulation Logic =================
    let r1_path = format!("{}/demo_R1.fastq", out_dir);
    let r2_path = format!("{}/demo_R2.fastq", out_dir);
    let mut r1_writer = BufWriter::new(File::create(r1_path)?);
    let mut r2_writer = BufWriter::new(File::create(r2_path)?);

    let num_pairs = (REF_LENGTH * COVERAGE) / (READ_LENGTH * 2);
    println!("Simulating {} million read pairs...", num_pairs / 1_000_000);

    for i in 0..num_pairs {
        // Randomly pick a starting position for the pair
        let p1_start = rng.gen_range(0..(REF_LENGTH - INSERT_SIZE - READ_LENGTH - 500));

        // Calculate the partner's position based on a normal distribution of insert size
        let current_insert = (INSERT_SIZE as i32 + rng.gen_range(-(INSERT_STD as i32)..INSERT_STD as i32)) as usize;
        let p2_start = p1_start + current_insert;

        let mut seq1 = reference[p1_start..p1_start + READ_LENGTH].to_vec();
        let mut seq2 = reference[p2_start..p2_start + READ_LENGTH].to_vec();

        // Introduce SVs into 1% of the read pairs
        let roll: f64 = rng.gen();
        if roll < 0.01 {
            // DELETION: Make the insert size appear huge to Lumpy
            let del_len = 1000;
            if p2_start + del_len + READ_LENGTH < REF_LENGTH {
                seq2 = reference[p2_start + del_len..p2_start + del_len + READ_LENGTH].to_vec();
            }
        } else if roll < 0.02 {
            // SOFT-CLIP: Cut the end of read 1 and replace with random DNA
            let clip_len = 50;
            for j in (READ_LENGTH - clip_len)..READ_LENGTH {
                seq1[j] = nucleotides[rng.gen_range(0..4)] as char;
            }
        }

        // Fastq Formatting
        let s1: String = seq1.iter().collect();
        let s2: String = seq2.iter().collect();
        // Reverse complement R2 to simulate Illumina Paired-End orientation (FR)
        let s2_rc: String = s2.chars().rev().map(|c| match c {
            'A'=>'T','T'=>'A','C'=>'G','G'=>'C',_=>'N'}).collect();

        writeln!(r1_writer, "@read{}/1\n{}\n+\n{}", i, s1, "I".repeat(READ_LENGTH))?;
        writeln!(r2_writer, "@read{}/2\n{}\n+\n{}", i, s2_rc, "I".repeat(READ_LENGTH))?;
    }

    println!("Simulation Complete. FASTQ files ready for alignment.");
    Ok(())
}