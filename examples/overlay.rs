//! `overlay <in.raw> <out.raw>` — real⊕sim: add sim peaks onto existing scans.
use thermorawfile::RawFile;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mut rf = RawFile::open(&a[1]).unwrap();
    let cal = rf.calibration_at_event(rf.scantrailer_addr as usize + 4).unwrap();

    // --- MS1 profile overlay (scan 1) ---
    let prof0 = rf.profile(1).unwrap();
    let (pts0, tic0) = (prof0.point_count(), prof0.chunks.iter().flat_map(|c| c.signal.iter()).map(|&v| v as f64).sum::<f64>());
    let sim_ms1 = [(555.5_f64, 9.9e5_f32), (744.4, 7.7e5)];
    rf.overlay_profile(1, &sim_ms1, &cal).unwrap();

    // --- MS2 centroid overlay (scan 2) ---
    let real2 = rf.centroid_peaks(2).len();
    let sim_ms2 = [(333.33_f64, 5.0e5_f32), (888.88, 4.0e5)];
    rf.overlay_centroids(2, &sim_ms2, 10.0).unwrap();

    rf.save(&a[2]).unwrap();

    // verify
    let rf2 = RawFile::open(&a[2]).unwrap();
    assert!(rf2.checksum_valid());
    let prof1 = rf2.profile(1).unwrap();
    let tic1: f64 = prof1.chunks.iter().flat_map(|c| c.signal.iter()).map(|&v| v as f64).sum();
    // sim m/z present?
    let has = |mz: f64| prof1.chunks.iter().any(|c| (0..c.signal.len()).any(|j| (prof1.mz_of_bin(c.first_bin + j as u32, &cal) - mz).abs() < 0.02));
    println!("MS1: real points {pts0} (TIC {tic0:.0}) -> overlaid points {} (TIC {tic1:.0})", prof1.point_count());
    println!("  555.5 present: {}, 744.4 present: {}  (TIC grew by {:.0}, expected ~{:.0})", has(555.5), has(744.4), tic1 - tic0, 9.9e5 + 7.7e5);
    let cents = rf2.centroid_peaks(2);
    println!("MS2: real centroids {real2} -> overlaid {}", cents.len());
}
