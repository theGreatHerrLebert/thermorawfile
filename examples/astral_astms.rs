//! `astral_astms <in.raw> <out.raw> [scan]` — prototype/verify the Astral ASTMS
//! 8-byte `{ f32 m/z, f32 intensity }` centroid write path.
//!
//! Authors a synthetic peak ramp into one ASTMS (MS2) scan, saves (recomputing
//! the scan-index stats + Adler-32), then re-opens and asserts the peaks read
//! back. A successful round-trip proves the 8-byte path: had the writer used the
//! 12-byte FTMS stride, it would have overrun the packet and the re-read would
//! mismatch. Run Thermo's RawFileReader on `out.raw` afterwards to cross-check.

use thermorawfile::{Peak, RawFile};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 3 {
        eprintln!("usage: astral_astms <in.raw> <out.raw> [scan=2]");
        std::process::exit(2);
    }
    let scan: u32 = a.get(3).and_then(|s| s.parse().ok()).unwrap_or(2);

    let mut rf = RawFile::open(&a[1]).expect("open input");
    let orig = rf.centroid_peaks(scan);
    let n = orig.len();
    assert!(n > 0, "scan {scan} has no centroid peaks (profile/MS1?)");
    println!(
        "scan {scan}: {n} existing centroid peaks, mz {:.4}..{:.4}",
        orig.first().unwrap().mz,
        orig.last().unwrap().mz
    );

    // Synthetic ascending ramp, SAME count (no offset rebuild needed).
    let ramp: Vec<Peak> = (0..n)
        .map(|i| Peak {
            mz: 200.0 + i as f64 * 0.25,
            intensity: 100.0 + i as f32 * 3.0,
        })
        .collect();
    rf.set_centroid_peaks(scan, &ramp).expect("rewrite ASTMS peaks");
    rf.save(&a[2]).expect("save output");
    println!("authored {n}-peak ramp into scan {scan} -> {}", a[2]);

    // Re-open the Rust-written file and verify.
    let rf2 = RawFile::open(&a[2]).expect("re-open output");
    assert!(rf2.checksum_valid(), "Adler-32 checksum invalid after save");
    let back = rf2.centroid_peaks(scan);
    assert_eq!(back.len(), n, "peak count changed on round-trip");

    // 8-byte form stores m/z as f32, so compare at f32 precision.
    let mut max_dmz = 0f64;
    for (got, want) in back.iter().zip(ramp.iter()) {
        max_dmz = max_dmz.max((got.mz - (want.mz as f32) as f64).abs());
        assert!((got.intensity - want.intensity).abs() < 1e-3, "intensity mismatch");
    }
    assert!(max_dmz < 1e-3, "m/z mismatch beyond f32 precision: {max_dmz}");

    println!(
        "round-trip OK: {n} peaks match (max |Δmz|={max_dmz:.2e}, f32 precision), checksum valid"
    );
    println!("first 3 read back: {:?}", &back[..3.min(back.len())]);
}
