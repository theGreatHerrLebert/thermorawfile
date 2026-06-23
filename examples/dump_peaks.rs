//! dump_peaks <in.raw> [iso_center] [iso_tol]
//! Prints, as CSV to stdout: scan, ms_order, iso_center, rt, mz, intensity
//! for MS2 scans whose isolation center is within iso_tol of iso_center
//! (default: all MS2 scans). Top 20 peaks per scan by intensity.
use thermorawfile::RawFile;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let rf = RawFile::open(&a[1]).expect("open");
    let target: Option<f64> = a.get(2).and_then(|s| s.parse().ok());
    let tol: f64 = a.get(3).and_then(|s| s.parse().ok()).unwrap_or(0.1);
    println!("scan,ms_order,iso_center,iso_width,ce,n_peaks,mz,intensity");
    for scan in 1..=rf.scan_count() as u32 {
        let ev = match rf.scan_event(scan) { Some(e) => e, None => continue };
        if ev.ms_order < 2 { continue; }
        if let Some(t) = target {
            if (ev.isolation_center - t).abs() > tol { continue; }
        }
        let mut peaks = rf.centroid_peaks(scan);
        if peaks.is_empty() { continue; }
        peaks.sort_by(|x,y| y.intensity.partial_cmp(&x.intensity).unwrap());
        let n = peaks.len();
        for p in peaks.iter().take(20) {
            println!("{scan},{},{:.4},{:.4},{:.2},{n},{:.5},{:.1}",
                ev.ms_order, ev.isolation_center, ev.isolation_width, ev.collision_energy, p.mz, p.intensity);
        }
    }
}
