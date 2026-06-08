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
