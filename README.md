# cluster_rs
🚧 Work in progress experimental read clustering tool suitable for Nanopore amplicon samples. Uses a sketch based minimizer approach along with Weighted Jaccard Index for calculating similarities between reads.

## Requirements
- Linux OS (Ubuntu 24.04.2)
- Rust >= 1.88.0

## Installation
Clone the repository or download the source code. Enter the cluster_rs directory and run:<br>
`cargo build --release`

The generated binary is available in `target/release/cluster_rs`.

## Usage
Run with:<br>
`cluster_rs --fastq <reads.fastq.gz> --outdir <outdir>`<br>


Optional arguments:
<pre>
<b>-s/--sketch_size</b> [200].

<b>-k/--kmer_size</b> [15]. Adjust according to error rate.

<b>-w/--window_size</b> [30]. Number of consecutive kmers to choose minimizer from.

<b>-j/--jaccard_distance</b> [0.1]. Min Jaccard distance to consider a read part of a cluster.
</pre>

## Read preprocessing
The algorithm is greedy, which means it is a good idea to filter and sort the reads prior (e.g., with [`fastq_rs`](https://github.com/OscarAspelin95/fastq_rs)).
