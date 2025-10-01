use clap::{Parser, Subcommand};
use log::*;
use rayon::prelude::*;

mod bins;
mod io_utils;
mod kmer;
mod suffix;
mod confirm;
mod mpi_wrap;
mod vcfout;
mod cluster;
//mod generate_demo_data;
mod merge_vcf;

#[derive(Parser)]
#[command(name = "svx", version, about = "Fast SNP/SV caller (Rust + MPI)")]
struct Cli {
    #[command(subcommand)]
    cmd: Command,
}

#[derive(Subcommand)]
enum Command {
    Index {
        #[arg(long)] ref_fa: String,
        #[arg(long, default_value_t = 150_000)] bin: usize,
        #[arg(long, default_value_t = 4_000)] pad: usize,
        #[arg(long)] out: String,
    },
    Call {
        #[arg(long)] ref_fa: String,
        #[arg(long)] bam: String,
        #[arg(long)] bins: String,
        #[arg(long, default_value = "out")] outdir: String,
        #[arg(long, default_value_t = 32)] band: usize,
        #[arg(long, default_value_t = 31)] k: usize,
        #[arg(long, default_value_t = 50)] w: usize,
        #[arg(long, default_value_t = false)] mpi: bool,
        #[arg(long)] threads: Option<usize>,
    },
    Merge {
        #[arg(long)] outdir: String,
        #[arg(long)] vcf: String,
        #[arg(long)] ref_fa: String,
    }
}

fn main() -> anyhow::Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.cmd {
        Command::Index { ref_fa, bin, pad, out } => {
            bins::build_bins(&ref_fa, bin, pad, &out)?;
            info!("bins written to {out}");
        }
        Command::Call { ref_fa, bam, bins: bins_path, outdir, band, k, w, mpi, threads } => {
            if let Some(t) = threads {
                rayon::ThreadPoolBuilder::new().num_threads(t).build_global().ok();
            }
            if mpi {
                mpi_wrap::run_mpi(&ref_fa, &bam, &bins_path, &outdir, band, k, w)?;
            } else {
                let bins = io_utils::read_bins(&bins_path)?;
                bins.into_par_iter().for_each(|b| {
                    if let Err(e) = process_bin(&ref_fa, &bam, &outdir, b, band, k, w) {
                        error!("bin failed: {e}");
                    }
                });
            }
        }
        Command::Merge { outdir, vcf, ref_fa } => {
            vcfout::merge_json_to_vcf(&outdir, &vcf, &ref_fa)?;
            info!("VCF written to {vcf}");
        }
    }
    Ok(())
}

pub fn process_bin(ref_fa: &str, bam_path: &str, outdir: &str, bin: bins::Bin,
                   band: usize, k: usize, w: usize) -> anyhow::Result<()> {
    use io_utils::*;
    use kmer::*;
    use suffix::*;
    use confirm::*;

    let ref_window = fetch_reference_window(ref_fa, &bin.chrom, bin.start, bin.end)?;
    let sfx = SuffixWorkspace::build(&ref_window);

    let mut reader = open_bam(bam_path)?;
    let reads = fetch_reads_overlapping(&mut reader, &bin.chrom, bin.start, bin.end)?;

    let mindex = MinimizerIndex::build(&ref_window, k, w);
    let hotspots = discover_hotspots(&reads, &mindex, &ref_window, k);

    let calls: Vec<_> = hotspots.into_par_iter()
        .filter_map(|hs| {
            match confirm_breakpoint(&bin, &ref_window, &sfx, &hs, band) {
                Ok(Some(call)) => Some(call),
                Ok(None) => None,
                Err(e) => { error!("confirm err: {e}"); None }
            }
        })
        .collect();

    vcfout::write_bin_json(outdir, &bin, &calls)?;
    Ok(())
}
