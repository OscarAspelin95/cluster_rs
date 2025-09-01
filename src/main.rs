use clap::Parser;

use simple_logger::SimpleLogger;

mod args;
use args::Args;

mod cluster;
use cluster::cluster;

#[macro_use]
mod utils;
use utils::Config;

fn main() {
    SimpleLogger::new()
        .with_level(log::LevelFilter::Info)
        .init()
        .unwrap();

    let args = Args::parse();

    let cfg = Config {
        sketch_size: args.sketch_size,
        kmer_size: args.kmer_size,
        window_size: args.window_size,
        jaccard_distance: args.jaccard_distance,
    };

    cluster(&args.fastq, &cfg, &args.outdir).unwrap();
}
