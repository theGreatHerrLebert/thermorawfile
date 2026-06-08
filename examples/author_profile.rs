//! `author_profile <in.raw> <out.raw>` — author a fully synthetic FTMS profile
//! with an arbitrary peak count at *exact* chosen m/z, then verify the round-trip.
//!
//! This is variable-count authoring: one single-bin chunk per peak, placed at the
//! exact grid bin via the inverted calibration, written within the scan's packet
//! budget with all offsets unchanged (no global rebuild). Cross-check with Thermo's
//! RawFileReader: it should report peaks at the chosen m/z with the chosen intensities.

use thermorawfile::RawFile;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 3 {
        eprintln!("usage: author_profile <in.raw> <out.raw>");
        std::process::exit(2);
    }
    let scan = 1u32;
    // Arbitrary count, arbitrary exact m/z (note: 5 peaks, different from the
    // original profile's thousands of points).
    let peaks: &[(f64, f32)] = &[
        (450.12345, 5.0e5),
        (500.0, 1.0e6),
        (618.314, 7.5e5),
        (777.777, 2.0e6),
        (880.5, 3.0e5),
    ];

    let mut rf = RawFile::open(&a[1]).expect("open input");
    let calib = rf
        .calibration_at_event(rf.scantrailer_addr as usize + 4)
        .expect("read calibration");
    let old = rf.profile(scan).unwrap().point_count();
    println!("calibration nparam={} b={:.4e} c={:.4e}", calib.nparam, calib.b, calib.c);
    println!("original profile points = {old}; authoring {} synthetic peaks", peaks.len());

    rf.author_profile(scan, peaks, &calib).expect("author profile");
    rf.save(&a[2]).expect("save");

    // Re-open and verify the authored peaks decode at the requested m/z.
    let rf2 = RawFile::open(&a[2]).expect("re-open");
    assert!(rf2.checksum_valid(), "checksum invalid");
    let cal2 = rf2.calibration_at_event(rf2.scantrailer_addr as usize + 4).unwrap();
    let prof2 = rf2.profile(scan).expect("profile gone");
    assert_eq!(prof2.chunks.len(), peaks.len(), "chunk count != peak count");
    println!("\nauthored {} chunks; read back:", prof2.chunks.len());
    let mut max_ppm = 0f64;
    for (ch, &(want_mz, want_int)) in prof2.chunks.iter().zip(peaks) {
        let got_mz = prof2.mz_of_bin(ch.first_bin, &cal2);
        let got_int = ch.signal[0];
        let ppm = (got_mz - want_mz).abs() / want_mz * 1e6;
        max_ppm = max_ppm.max(ppm);
        println!(
            "  want m/z {want_mz:.5} int {want_int:.0}  ->  bin {} = m/z {got_mz:.5} int {got_int:.0}  ({ppm:.2} ppm)",
            ch.first_bin
        );
        assert!((got_int - want_int).abs() < 1.0, "intensity mismatch");
    }
    println!("\nround-trip OK: {} peaks, max {max_ppm:.2} ppm (grid quantization), checksum valid", peaks.len());
}
