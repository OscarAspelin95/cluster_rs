use anyhow::Result;
use fxhash::FxHasher;
use log::info;
use minimizer_iter::MinimizerBuilder;
use minimizer_iter::iterator::MinimizerIterator;
use needletail::parse_fastx_file;
use probminhash::jaccard::compute_probminhash_jaccard;
use probminhash::probminhasher::ProbMinHash2;
use std::collections::HashSet;
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use crate::Config;
use crate::file_path;

/// Return a lazy minimizer iterator.
fn get_minimizers(seq: &[u8], kmer_size: usize, window_size: u16) -> MinimizerIterator<'_> {
    let builder: MinimizerBuilder<u64> = MinimizerBuilder::new()
        .minimizer_size(kmer_size)
        .width(window_size);

    let minimizer_iterator: minimizer_iter::iterator::MinimizerIterator<'_> = builder.iter(seq);

    minimizer_iterator
}

/// Consider how minimizer counts should affect the weights.
/// * Low minimizer duplicity means rare.
/// * High minimizer duplicity means common.
fn sketch_read(seq: &[u8], cfg: &Config) -> Vec<u64> {
    let mut pmh = ProbMinHash2::<u64, FxHasher>::new(cfg.sketch_size, u64::MAX);

    let mut hmap: HashSet<u64> =
        HashSet::with_capacity(seq.len() - cfg.kmer_size - cfg.window_size as usize + 2);

    for (hash, _) in get_minimizers(seq, cfg.kmer_size, cfg.window_size) {
        if hmap.contains(&hash) {
            continue;
        }
        pmh.hash_item(hash, 1.0);
        hmap.insert(hash);
    }

    pmh.get_signature().to_owned()
}

/// * Consider replacing probminhash crate with manual implementation
/// to regain more control over clustering.
///
/// * Consider testing a global/semi-global aligner instead of minimizers.
pub fn cluster(fastq: &PathBuf, cfg: &Config, outdir: &PathBuf) -> Result<()> {
    if !outdir.is_dir() {
        create_dir_all(&outdir)?;
    }

    let mut reader = parse_fastx_file(&fastq)?;

    info!("Sketching reads...");

    // cluster_id, cluster_sketch, num_reads_belonging.
    let mut clusters: Vec<(usize, Vec<u64>, usize)> = Vec::new();

    while let Some(record) = reader.next() {
        let record = match record {
            Ok(record) => record,
            Err(_) => continue,
        };

        let seq = record.seq();

        // Minimum length to accomodate w consecutive kmers.
        let required_length: usize = cfg.window_size as usize + cfg.kmer_size - 1;
        if seq.len() < required_length {
            continue;
        }

        // MinHash sketch for read.
        let read_sketch = sketch_read(&seq, &cfg);

        // Check if read belongs to an already existing cluster.
        let mut assigned: bool = false;

        for (_, cluster_sketch, cluster_members) in &mut clusters {
            let jaccard_distance =
                compute_probminhash_jaccard(&read_sketch[..], &cluster_sketch[..]);

            // Read belongs to an existing cluster.
            if jaccard_distance >= cfg.jaccard_distance {
                *cluster_members += 1;
                assigned = true;
                break;
            }
        }

        // Initialize new cluster.
        if !assigned {
            clusters.push((clusters.len() + 1, read_sketch, 1));

            let representative_file =
                file_path!(&outdir, format!("cluster_{}.fastq", clusters.len()));

            let mut writer = BufWriter::new(File::create(representative_file)?);

            // Id.
            writer.write_all(b">")?;
            writer.write_all(record.id())?;

            // Sequence.
            writer.write_all(b"\n")?;
            writer.write_all(&record.seq())?;
        }
    }

    // Sort clusters by size.
    clusters.sort_by(|a, b| a.2.cmp(&b.2).reverse());

    // Write cluster info to file.
    let cluster_file = file_path!(outdir, "clusters.tsv");
    let mut writer = BufWriter::new(File::create(cluster_file)?);

    writer.write_all(b"cluster_id\tcluster_members\n")?;

    for (cid, _, cmembers) in clusters {
        writer.write_all(cid.to_string().as_bytes())?;
        writer.write_all(b"\t")?;
        writer.write_all(cmembers.to_string().as_bytes())?;
        writer.write_all(b"\n")?;
    }

    Ok(())
}
