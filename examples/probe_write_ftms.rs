//! Round-trip authoring into a real Orbitrap (FTMS 12-byte centroid) MS2 slot:
//! prove the writer emits the analyzer-correct centroid format + valid checksum.
use thermorawfile::RawFile;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mut rf = RawFile::open(&a[1]).expect("open");
    // first MS2 scan
    let mut ms2 = 0u32;
    for s in 1..=rf.scan_count() as u32 {
        if let Some(e) = rf.scan_event(s) { if e.ms_order >= 2 { ms2 = s; break; } }
    }
    let before = rf.centroid_peaks(ms2).len();
    // author 3 known peaks (high-precision m/z to exercise the f64 wide path)
    let peaks = [(345.67891_f64, 1000.0_f32), (678.12345, 500.0), (912.98765, 250.0)];
    rf.author_centroids(ms2, &peaks).expect("author");
    rf.save(&a[2]).expect("save");
    let rf2 = RawFile::open(&a[2]).expect("reopen");
    let got = rf2.centroid_peaks(ms2);
    println!("MS2 scan {ms2}: before={before} peaks -> authored 3, read back {}", got.len());
    println!("checksum_valid={}", rf2.checksum_valid());
    for (i, p) in got.iter().take(3).enumerate() {
        println!("  peak{i}: mz={:.5} int={:.1}  (target {:.5})", p.mz, p.intensity, peaks[i].0);
    }
    // m/z precision check: wide (FTMS) stores f64; narrow (ASTMS) would truncate to f32
    let err = (got[0].mz - peaks[0].0).abs();
    println!("first-peak m/z abs error = {:.2e}  ({})", err,
        if err < 1e-6 { "f64/FTMS wide OK" } else { "f32/narrow (truncated)" });
}
