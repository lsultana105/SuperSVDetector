# SuperSVDetector

SuperSVDetector is a Rust-based structural variant (SV) detection pipeline designed to improve the computational efficiency of short-read whole-genome analysis. The system focuses on scalable detection of structural variants, particularly deletions, using a bin-based processing strategy combined with both shared-memory and distributed (MPI) parallelism.

This project was developed as part of a Master's dissertation investigating whether a Rust-based implementation can reproduce validated SVs while significantly improving runtime performance compared to existing tools such as LUMPY.

---

## Features

- High-performance implementation in Rust
- Genome partitioning into independent bins
- Hotspot discovery using minimizer-based k-mer indexing
- Breakpoint confirmation using alignment evidence
- MPI-based distributed execution (master–worker model)
- Per-bin JSON output with final VCF merging
- Designed for large-scale BAM file processing

---

## Pipeline Overview

The SuperSVDetector pipeline consists of the following stages:

1. Partition the reference genome into bins
2. Process each bin independently
3. Identify candidate hotspots using read evidence
4. Confirm breakpoints using local alignment features
5. Output per-bin structural variant calls in JSON format
6. Merge all results into a final VCF file

---

## Repository Structure
src/
├── main.rs
├── process_bin.rs
├── bins.rs
├── io_utils.rs
├── kmer.rs
├── hotspot.rs
├── confirm.rs
├── mpi_wrap.rs
├── cluster.rs
└── vcfout.rs


---

## Requirements

- Rust (latest stable)
- Clang (required for building dependencies)
- OpenMPI (for distributed execution)
- samtools (for indexing BAM and FASTA files)

Input requirements:
- BAM file with `.bai` index
- Reference FASTA file with `.fai` index
- Matching chromosome naming (e.g., `chr1` vs `1`)

---

## Build Instructions

```bash
export CC=clang
export CXX=clang++
cargo clean
cargo build --release

Usage
**Step 1: Generate bins
**For example:
./target/release/SuperSVDetector index \
  --ref-fa reference.fa \
  --bin-size 1000000 \
  --out bins.tsv

**Step 2: Run variant detection**
./target/release/SuperSVDetector call \
  --ref-fa reference.fa \
  --bam sample.bam \
  --bins bins.tsv \
  --outdir out_calls \
  --k 31 \
  --w 32 \
  --band 500

**Step 3: Merge results into VCF**
./target/release/SuperSVDetector merge \
  --json-dir out_calls \
  --out variants.vcf

MPI Execution
mpirun -n 4 ./target/release/SuperSVDetector call \
  --mpi \
  --ref-fa reference.fa \
  --bam sample.bam \
  --bins bins.tsv \
  --outdir out_calls

