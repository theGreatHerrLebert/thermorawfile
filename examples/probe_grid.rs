use thermorawfile::RawFile;
fn main() {
    let rf = RawFile::open(format!("{}/tests/data/small2.RAW", env!("CARGO_MANIFEST_DIR"))).unwrap();
    let cal = rf.calibration_at_event(rf.scantrailer_addr as usize + 4).unwrap();
    let prof = rf.profile(1).unwrap();
    let (fv, step, nbins) = (prof.first_value, prof.step, prof.nbins);
    println!("fv={fv} step={step} nbins={nbins} point_count={}", prof.point_count());
    println!("centroid_peaks(1)={}", rf.centroid_peaks(1).len());
    for bin in [0u32, 1, 100, nbins / 4, nbins / 2, nbins - 2] {
        let mz = cal.mz(fv + bin as f64 * step);
        let back = cal.freq(mz).map(|f| ((f - fv) / step).round());
        println!("bin {bin:>8} -> mz {mz:.4} -> recovered_bin {back:?}");
    }
}
