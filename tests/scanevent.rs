//! Scan-event (acquisition descriptor) read + author round-trip, cross-checked
//! against values from Thermo's RawFileReader.

use thermorawfile::RawFile;

fn sample(name: &str) -> String {
    format!("{}/tests/data/{}", env!("CARGO_MANIFEST_DIR"), name)
}

#[test]
fn reads_scan_events() {
    let rf = RawFile::open(sample("small2.RAW")).unwrap();
    assert_eq!(rf.scan_event_size, 232, "rev-66 fixed event stride");
    // scan 1 = MS1 (FTMS), scan 2 = MS2 (ITMS) with a precursor.
    let e1 = rf.scan_event(1).unwrap();
    assert_eq!(e1.ms_order, 1);
    assert_eq!(e1.analyzer, 4); // FTMS
    let e2 = rf.scan_event(2).unwrap();
    assert_eq!(e2.ms_order, 2);
    assert_eq!(e2.analyzer, 0); // ITMS
    // Oracle: scan 2 precursor 398.5411, CE 35.0.
    assert!((e2.isolation_center - 398.5411).abs() < 1e-3, "got {}", e2.isolation_center);
    assert!((e2.collision_energy - 35.0).abs() < 1e-6);
}

#[test]
fn authors_isolation_window() {
    let mut rf = RawFile::open(sample("small2.RAW")).unwrap();
    rf.set_isolation(2, 555.5, 13.0, 27.5).unwrap();
    let out = format!("{}/event.raw", env!("CARGO_TARGET_TMPDIR"));
    rf.save(&out).unwrap();

    let rf2 = RawFile::open(&out).unwrap();
    assert!(rf2.checksum_valid());
    let e = rf2.scan_event(2).unwrap();
    assert!((e.isolation_center - 555.5).abs() < 1e-9);
    assert!((e.isolation_width - 13.0).abs() < 1e-9);
    assert!((e.collision_energy - 27.5).abs() < 1e-9);
    assert_eq!(e.ms_order, 2); // unchanged
    assert_eq!(e.analyzer, 0); // unchanged
}
