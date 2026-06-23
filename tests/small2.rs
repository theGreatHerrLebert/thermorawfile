//! Validates the pure-Rust reader against the real rev-66 Orbitrap sample from
//! the public ThermoRawFileParser corpus. Expected peaks are reference values
//! present in that public file (self-consistent with each scan's own structure).

use thermorawfile::{compute_checksum, stored_checksum, RawFile};

fn sample(name: &str) -> String {
    format!("{}/tests/data/{}", env!("CARGO_MANIFEST_DIR"), name)
}

#[test]
fn reads_rev66_structure() {
    let rf = RawFile::open(sample("small2.RAW")).expect("open small2.RAW");
    assert_eq!(rf.version, 66);
    assert_eq!(rf.first_scan, 1);
    assert_eq!(rf.last_scan, 95);
    assert_eq!(rf.scan_count(), 95);
}

#[test]
fn checksum_matches_both_revisions() {
    for f in ["small2.RAW", "small.RAW"] {
        let bytes = std::fs::read(sample(f)).expect("read sample");
        assert_eq!(
            compute_checksum(&bytes),
            stored_checksum(&bytes),
            "checksum mismatch for {f}"
        );
    }
}

#[test]
fn reads_scan2_centroid_peaks() {
    let rf = RawFile::open(sample("small2.RAW")).expect("open");
    let peaks = rf.centroid_peaks(2);
    assert_eq!(peaks.len(), 196, "scan 2 peak count");
    // Reference peaks from the public small2.RAW fixture: [0] m/z=116.0264 int=12.1, [195] m/z=882.6021.
    assert!((peaks[0].mz - 116.0264).abs() < 1e-3, "got {}", peaks[0].mz);
    assert!((peaks[0].intensity - 12.1).abs() < 0.5);
    assert!((peaks[195].mz - 882.6021).abs() < 1e-3, "got {}", peaks[195].mz);
}

#[test]
fn reads_scan2_trailer_params() {
    let rf = RawFile::open(sample("small2.RAW")).expect("open");
    let p = rf.scan_params(2).expect("scan 2 trailer params");
    assert_eq!(p.charge_state(), Some(3));
    assert_eq!(p.ion_injection_time_ms(), Some(50.0));
    // cross-check: the trailer's monoisotopic m/z equals the scan-2 precursor center.
    assert!((p.monoisotopic_mz().unwrap() - 398.5411).abs() < 1e-3);
    assert_eq!(p.isolation_width_mz(), Some(2.0));
    assert_eq!(p.micro_scan_count(), Some(1));
    // the underlying record exposes every label, not just the typed accessors
    assert!(p.record().get("FT Resolution:").is_some());
}
