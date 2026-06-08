//! `author_ms2 <in.raw> <out.raw>` — author arbitrary-count ASTMS MS2 centroids.
use thermorawfile::RawFile;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let scan = 2u32; // first ASTMS MS2 scan in the Astral sample
    let peaks: &[(f64, f32)] = &[(120.0805, 4.0e4), (175.119, 9.0e4), (550.27, 2.5e5), (876.5, 1.1e5)];
    let mut rf = RawFile::open(&a[1]).unwrap();
    rf.author_centroids(scan, peaks).expect("author centroids");
    rf.save(&a[2]).unwrap();
    let rf2 = RawFile::open(&a[2]).unwrap();
    assert!(rf2.checksum_valid());
    let back = rf2.centroid_peaks(scan);
    println!("authored {} MS2 peaks; read back {}:", peaks.len(), back.len());
    for p in &back { println!("  m/z={:.4} int={:.0}", p.mz, p.intensity); }
}
