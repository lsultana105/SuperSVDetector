use std::fs::File;
use std::io::{Write, BufWriter};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

// ==== Constants for Accuracy Testing ====
const REF_LENGTH: usize = 5_000_000;   // 5 Mb toy ref
const LINE_WIDTH: usize = 80;
const READ_LENGTH: usize = 150;
const COVERAGE: usize = 10;          // also lower this initially
const BIN_SIZE: usize = 1_000_000;     // 1 Mb bins
const INSERT_SIZE: usize = 450;        // Distance between paired reads
const INSERT_STD: usize = 20;

fn main() -> std::io::Result<()> {
    let out_dir = "dataToTest";
    std::fs::create_dir_all(out_dir)?;

    // ================= 1. Define Ground Truth =================
    // 8 Deletions and 5 Insertions at specific points
    // Put deletions and insertions clearly inside this range
    let deletions: Vec<usize> = vec![
        500_000,   // 0.5 Mb
        1_000_000, // 1.0 Mb
        1_500_000, // 1.5 Mb
        2_000_000, // 2.0 Mb
    ];

    let insertions: Vec<usize> = vec![
        750_000,   // between 0.5 and 1.0 Mb
        1_250_000, // between 1.0 and 1.5 Mb
        1_750_000, // between 1.5 and 2.0 Mb
    ];

    let deletion_len: usize = 2_000;
    let insertion_len: usize = 50;

    // ================= 2. Generate Reference =================
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
    println!("Reference (50MB) written to {}", ref_path);

    // ================= 3. Generate Bins =================
    let bins_path = format!("{}/demo_bins.txt", out_dir);
    let mut bins_file = BufWriter::new(File::create(bins_path)?);
    let mut start = 0;
    while start < REF_LENGTH {
        let end = std::cmp::min(start + BIN_SIZE, REF_LENGTH);
        writeln!(bins_file, "chr_sim_1\t{}\t{}", start, end)?;
        start = end;
    }

    // ================= 4. Simulation Logic =================
    let r1_path = format!("{}/demo_R1.fastq", out_dir);
    let r2_path = format!("{}/demo_R2.fastq", out_dir);
    let mut r1_writer = BufWriter::new(File::create(r1_path)?);
    let mut r2_writer = BufWriter::new(File::create(r2_path)?);

    let num_pairs = (REF_LENGTH * COVERAGE) / (READ_LENGTH * 2);
    println!("Simulating {} read pairs for Accuracy Test...", num_pairs);

    for i in 0..num_pairs {
        // Pick a random starting position
        let p1_start = rng.gen_range(500..(REF_LENGTH - 2000));
        let mut p2_offset = (INSERT_SIZE as i32 + rng.gen_range(-(INSERT_STD as i32)..INSERT_STD as i32)) as usize;

        let mut seq1 = reference[p1_start..p1_start + READ_LENGTH].to_vec();

        // --- INJECT DELETION SIGNAL ---
        // If the read pair straddles a deletion site, push R2 further away
        for &pos in &deletions {
            if p1_start < pos && (p1_start + p2_offset) > pos {
                p2_offset += 2000; // Simulate a 2kb deletion
            }
        }

        let p2_start = p1_start + p2_offset;
        let mut seq2 = reference[p2_start..p2_start + READ_LENGTH].to_vec();

        // --- INJECT INSERTION (SOFT-CLIP) SIGNAL ---
        // If a read starts exactly at an insertion site, clip it
        for &pos in &insertions {
            if p1_start > pos - 10 && p1_start < pos + 10 {
                // Replace last 50bp with "Foreign DNA" (all A's)
                for j in (READ_LENGTH - 50)..READ_LENGTH {
                    seq1[j] = 'A';
                }
            }
        }

        // Fastq Formatting
        let s1: String = seq1.iter().collect();
        let s2: String = seq2.iter().collect();
        let s2_rc: String = s2.chars().rev().map(|c| match c {
            'A'=>'T','T'=>'A','C'=>'G','G'=>'C',_=>'N'}).collect();

        writeln!(r1_writer, "@read{}/1\n{}\n+\n{}", i, s1, "I".repeat(READ_LENGTH))?;
        writeln!(r2_writer, "@read{}/2\n{}\n+\n{}", i, s2_rc, "I".repeat(READ_LENGTH))?;
    }

    println!("Simulation Complete. Truth set: 8 Deletions, 5 Insertions.");
    Ok(())
}