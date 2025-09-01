use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct Args {
    #[arg(short, long, help = "Path to fastq file.")]
    pub fastq: PathBuf,

    #[arg(short, long, default_value_t = 200)]
    pub sketch_size: usize,

    #[arg(short, long, default_value_t = 15)]
    pub kmer_size: usize,

    #[arg(short, long, default_value_t = 30)]
    pub window_size: u16,

    #[arg(short, long, default_value_t = 0.1)]
    pub jaccard_distance: f64,

    #[arg(short, long, help = "Outdir.")]
    pub outdir: PathBuf,
}
