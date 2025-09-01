#[macro_export]
macro_rules! file_path {
    ($base:expr $(, $sub:expr)+) => {{
        use ::std::path::PathBuf;

        let mut full_path = PathBuf::from($base);

        $(
            full_path.push($sub);
        )*

        full_path
    }};
}

pub struct Config {
    pub sketch_size: usize,
    pub kmer_size: usize,
    pub window_size: u16,
    pub jaccard_distance: f64,
}
