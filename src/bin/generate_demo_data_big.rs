use std::fs::File;
use std::io::{Write, BufWriter};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

// ==== Larger Synthetic Genome for Runtime Testing ====
// Aim: stress-test runtime & MPI, not perfect biology.
const REF_LENGTH: usize = 200_000_000; // 200 Mb (4x your 50 Mb test)
const LINE_WIDTH: usize = 80;
const READ_LENGTH: usize = 150;
const COVERAGE: usize = 30;           // same depth as before
const BIN_SIZE: usize = 1_000_000;    // 1 Mb bins
const INSERT_SIZE: usize = 450;       // same as before
const INSERT_STD: usize = 20;

fn main() -> std::io::Result<()> {
    let out_dir = "dataToTestBig";
    std::fs::create_dir_all(out_dir)?;

    // ================= 1. Define Ground Truth =================
    // Place deletions every ~20 Mb, insertions in between.
    let mut deletions = Vec::new();
    let mut insertions = Vec::new();

    let step = 20_000_000; // 20 Mb spacing
    let del_size = 2000usize;
    let ins_len = 50usize;

    let mut pos = 10_000_000; // start at 10 Mb to avoid edges
    while pos + step < REF_LENGTH {
        deletions.push(pos);
        // insertion halfway between deletions
        insertions.push(pos + step / 2);
        pos += step;
    }

    println!("Planned {} deletions and {} insertions.", deletions.len(), insertions.len());

    // ================= 2. Generate Reference =================
    let ref_path = format!("{}/big_ref.fa", out_dir);
    let mut ref_writer = BufWriter::new(File::create(&ref_path)?);
    writeln!(ref_writer, ">chr_big_1")?;

    let mut rng = StdRng::seed_from_u64(42);
    let nucleotides = b"ACGT";
    let mut reference: Vec<char> = Vec::with_capacity(REF_LENGTH);

    for _ in 0..REF_LENGTH {
        reference.push(nucleotides[rng.gen_range(0..4)] as char);
    }

    for (i, &nt) in reference.iter().enumerate() {
        write!(ref_writer, "{}", nt)?;
        if (i + 1) % LINE_WIDTH == 0 {
            writeln!(ref_writer)?;
        }
    }
    writeln!(ref_writer)?;
    println!("Reference (200Mb) written to {}", ref_path);

    // ================= 3. Generate Bins =================
    let bins_path = format!("{}/big_bins.txt", out_dir);
    let mut bins_file = BufWriter::new(File::create(&bins_path)?);
    let mut start = 0;
    while start < REF_LENGTH {
        let end = std::cmp::min(start + BIN_SIZE, REF_LENGTH);
        writeln!(bins_file, "chr_big_1\t{}\t{}", start, end)?;
        start = end;
    }
    println!("Bin file written to {}", bins_path);

    // ================= 4. FASTQ Simulation =================
    let r1_path = format!("{}/big_R1.fastq", out_dir);
    let r2_path = format!("{}/big_R2.fastq", out_dir);
    let mut r1_writer = BufWriter::new(File::create(r1_path)?);
    let mut r2_writer = BufWriter::new(File::create(r2_path)?);

    let num_pairs = (REF_LENGTH * COVERAGE) / (READ_LENGTH * 2);
    println!("Simulating {} read pairs for runtime test...", num_pairs);

    for i in 0..num_pairs {
        // Pick a random starting position with some safety margin for insert size.
        let p1_start = rng.gen_range(1000..(REF_LENGTH - 5000));
        let mut p2_offset = (INSERT_SIZE as i32
            + rng.gen_range(-(INSERT_STD as i32)..(INSERT_STD as i32))) as usize;

        let mut seq1 = reference[p1_start..p1_start + READ_LENGTH].to_vec();

        // --- DELETION SIGNAL: push mate further away if spanning a deletion ---
        for &pos in &deletions {
            if p1_start < pos && (p1_start + p2_offset) > pos {
                p2_offset += del_size; // emulate ~2kb deletion
            }
        }

        let p2_start = p1_start + p2_offset;
        if p2_start + READ_LENGTH >= REF_LENGTH {
            continue;
        }

        let mut seq2 = reference[p2_start..p2_start + READ_LENGTH].to_vec();

        // --- INSERTION SIGNAL: soft-clip reads around insertion sites ---
        for &pos in &insertions {
            if p1_start > pos - 10 && p1_start < pos + 10 {
                // clip last 50bp to 'A's like before
                for j in (READ_LENGTH - ins_len)..READ_LENGTH {
                    seq1[j] = 'A';
                }
            }
        }

        let s1: String = seq1.iter().collect();
        let s2: String = seq2.iter().collect();
        let s2_rc: String = s2.chars().rev().map(|c| match c {
            'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N'
        }).collect();

        writeln!(r1_writer, "@read{}/1\n{}\n+\n{}", i, s1, "I".repeat(READ_LENGTH))?;
        writeln!(r2_writer, "@read{}/2\n{}\n+\n{}", i, s2_rc, "I".repeat(READ_LENGTH))?;
    }

    println!("Big simulation complete. Deletions: {}, Insertions: {}.", deletions.len(), insertions.len());
    Ok(())
}
