//! Pure-Rust write round-trip: read -> rewrite a scan's peaks -> save -> re-read,
//! and assert the checksum validates and the peaks come back exactly.

use thermorawfile::{Peak, RawFile};

fn sample(name: &str) -> String {
    format!("{}/tests/data/{}", env!("CARGO_MANIFEST_DIR"), name)
}

fn synthetic(n: usize) -> Vec<Peak> {
    (0..n)
        .map(|i| Peak {
            mz: 250.0 + i as f64 * 0.5,
            intensity: 700.0 + i as f32 * 4.0,
        })
        .collect()
}

fn synthetic_pairs(n: usize) -> Vec<(f64, f32)> {
    synthetic(n).into_iter().map(|p| (p.mz, p.intensity)).collect()
}

#[test]
fn pure_rust_write_roundtrip() {
    let mut rf = RawFile::open(sample("small2.RAW")).unwrap();
    let n = rf.centroid_peaks(2).len();
    assert_eq!(n, 196);

    let peaks = synthetic(n);
    rf.set_centroid_peaks(2, &peaks).unwrap();

    let out = format!("{}/written.raw", env!("CARGO_TARGET_TMPDIR"));
    rf.save(&out).unwrap();

    // Re-open the file we just wrote, with our own reader.
    let rf2 = RawFile::open(&out).unwrap();
    assert!(rf2.checksum_valid(), "Adler-32 must validate after write");

    let got = rf2.centroid_peaks(2);
    assert_eq!(got.len(), n);
    assert!((got[0].mz - 250.0).abs() < 1e-6, "got {}", got[0].mz);
    assert!((got[0].intensity - 700.0).abs() < 1e-3);
    assert!((got[n - 1].mz - (250.0 + (n - 1) as f64 * 0.5)).abs() < 1e-6);

    // A different scan must be untouched.
    let orig3 = RawFile::open(sample("small2.RAW")).unwrap().centroid_peaks(3);
    let new3 = rf2.centroid_peaks(3);
    assert_eq!(orig3.len(), new3.len());
    assert!((orig3[0].mz - new3[0].mz).abs() < 1e-9);
}

#[test]
fn rejects_count_mismatch() {
    let mut rf = RawFile::open(sample("small2.RAW")).unwrap();
    // wrong count -> error (variable-count rewrite is TODO)
    assert!(rf.set_centroid_peaks(2, &synthetic(10)).is_err());
}

/// Tier 1: repack a centroid scan to MORE peaks than its original packet budget,
/// then re-read with our own reader and assert the whole file stays consistent.
#[test]
fn repack_grow_roundtrip() {
    let mut rf = RawFile::open(sample("small2.RAW")).unwrap();
    let orig = rf.centroid_peaks(2).len();
    assert_eq!(orig, 196);

    // 3x the original peak count — author_centroids would reject this as over-budget.
    let grown = synthetic_pairs(orig * 3);
    assert!(
        rf.author_centroids(2, &grown).is_err(),
        "precondition: in-place author must reject the over-budget write"
    );
    rf.repack_centroids(2, &grown).unwrap();

    let out = format!("{}/repacked_grow.raw", env!("CARGO_TARGET_TMPDIR"));
    rf.save(&out).unwrap();

    let rf2 = RawFile::open(&out).unwrap();
    assert!(rf2.checksum_valid(), "Adler-32 must validate after repack");

    // The grown scan reads back exactly.
    let got = rf2.centroid_peaks(2);
    assert_eq!(got.len(), grown.len());
    assert!((got[0].mz - 250.0).abs() < 1e-6);
    assert!((got[got.len() - 1].mz - (250.0 + (grown.len() - 1) as f64 * 0.5)).abs() < 1e-6);

    // Every OTHER scan is byte-for-byte intact (peaks unchanged after the relocation).
    let base = RawFile::open(sample("small2.RAW")).unwrap();
    for s in base.first_scan..=base.last_scan {
        if s == 2 {
            continue;
        }
        let a = base.centroid_peaks(s);
        let b = rf2.centroid_peaks(s);
        assert_eq!(a.len(), b.len(), "scan {s} peak count changed");
        for (x, y) in a.iter().zip(b.iter()) {
            assert!((x.mz - y.mz).abs() < 1e-9 && (x.intensity - y.intensity).abs() < 1e-3,
                "scan {s} peak drift after repack");
        }
    }

    // Scan filters (built from the relocated scan-event stream) still parse.
    assert!(rf2.scan_filter(base.first_scan).is_some());
    assert!(rf2.scan_filter(base.last_scan).is_some());
}

/// Repack to FEWER peaks (shrink) must also keep the file consistent.
#[test]
fn repack_shrink_roundtrip() {
    let mut rf = RawFile::open(sample("small2.RAW")).unwrap();
    let small = synthetic_pairs(5);
    rf.repack_centroids(2, &small).unwrap();
    let out = format!("{}/repacked_shrink.raw", env!("CARGO_TARGET_TMPDIR"));
    rf.save(&out).unwrap();
    let rf2 = RawFile::open(&out).unwrap();
    assert!(rf2.checksum_valid());
    assert_eq!(rf2.centroid_peaks(2).len(), 5);
    // A later scan still intact.
    let base = RawFile::open(sample("small2.RAW")).unwrap();
    assert_eq!(base.centroid_peaks(3).len(), rf2.centroid_peaks(3).len());
}
