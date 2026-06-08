//! `ftms_profile <in.raw> <out.raw> [scan]` — prototype/verify the MS1 FTMS
//! profile encode path.
//!
//! Decodes an FTMS profile (chunked frequency-grid signal), overwrites every
//! signal point with a known value in place (grid preserved), saves (recomputing
//! TIC + Adler-32), then re-opens and asserts the profile reads back. This is
//! the "rewrite intensities on a real m/z grid" primitive for emitting a
//! simulated MS1 onto an Astral/Orbitrap template. Cross-check with Thermo's
//! RawFileReader: the scan's reported TIC must equal point_count * fill value.

use thermorawfile::RawFile;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 3 {
        eprintln!("usage: ftms_profile <in.raw> <out.raw> [scan=1]");
        std::process::exit(2);
    }
    let scan: u32 = a.get(3).and_then(|s| s.parse().ok()).unwrap_or(1);
    const FILL: f32 = 1000.0;

    let mut rf = RawFile::open(&a[1]).expect("open input");
    let prof = rf.profile(scan).expect("scan has no FTMS profile (centroid-only?)");
    let n = prof.point_count();
    println!(
        "scan {scan}: FTMS profile — {} chunks, {n} signal points, first_value={:.4} step={:.3e}, nbins={}",
        prof.chunks.len(),
        prof.first_value,
        prof.step,
        prof.nbins
    );

    // Overwrite every signal point with a constant; expected TIC = n * FILL.
    let new_signal = vec![FILL; n];
    rf.set_profile_intensities(scan, &new_signal)
        .expect("rewrite profile intensities");
    rf.save(&a[2]).expect("save output");
    println!("filled {n} profile points with {FILL} into scan {scan} -> {}", a[2]);
    println!("expected scan TIC = {n} * {FILL} = {:.0}", n as f64 * FILL as f64);

    // Re-open and verify.
    let rf2 = RawFile::open(&a[2]).expect("re-open output");
    assert!(rf2.checksum_valid(), "Adler-32 invalid after save");
    let prof2 = rf2.profile(scan).expect("profile gone after rewrite");
    assert_eq!(prof2.point_count(), n, "profile point count changed");
    assert_eq!(prof2.chunks.len(), prof.chunks.len(), "chunk count changed");
    let all_match = prof2
        .chunks
        .iter()
        .flat_map(|c| c.signal.iter())
        .all(|&v| (v - FILL).abs() < 1e-3);
    assert!(all_match, "profile signal did not read back as the fill value");
    // grid must be byte-preserved
    assert_eq!(prof2.chunks[0].first_bin, prof.chunks[0].first_bin, "grid changed");
    assert!((prof2.step - prof.step).abs() < 1e-12, "step changed");

    println!("round-trip OK: {n} profile points read back as {FILL}, grid preserved, checksum valid");
}
