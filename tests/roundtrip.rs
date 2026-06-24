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

/// Tier 1 (profile): repack an FTMS MS1 profile scan to MORE peaks than its budget,
/// binning onto the existing grid, then re-read and assert consistency.
#[test]
fn repack_profile_grow_roundtrip() {
    let mut rf = RawFile::open(sample("small2.RAW")).unwrap();
    // Scan 1 is the MS1 profile scan; its calibration is at scantrailer_addr + 4.
    let cal = rf
        .calibration_at_event(rf.scantrailer_addr as usize + 4)
        .expect("MS1 calibration");
    let prof = rf.profile(1).expect("scan 1 has a profile");
    let (fv, step, nbins) = (prof.first_value, prof.step, prof.nbins);

    // Build K on-grid peaks, keeping only bins that round-trip to themselves through
    // the quadratic calibration (its freq() inverse is double-valued for some m/z —
    // a calibration property, independent of the repack). K > the scan's 3032 profile
    // points, so the packet must grow.
    let k = 4000usize;
    let mut grown: Vec<(f64, f32)> = Vec::with_capacity(k);
    let mut bin = 1u32;
    while grown.len() < k && bin < nbins {
        let mz = cal.mz(fv + bin as f64 * step);
        if let Some(f) = cal.freq(mz) {
            let rb = ((f - fv) / step).round();
            if rb >= 0.0 && rb < nbins as f64 && rb as u32 == bin {
                let n = grown.len();
                grown.push((mz, 1000.0 + n as f32));
            }
        }
        bin += 1;
    }
    assert!(grown.len() >= 3500, "expected a grow; only {} valid bins", grown.len());
    let k = grown.len();

    rf.repack_profile(1, &grown, &cal).unwrap();
    let out = format!("{}/repacked_profile.raw", env!("CARGO_TARGET_TMPDIR"));
    rf.save(&out).unwrap();

    let rf2 = RawFile::open(&out).unwrap();
    assert!(rf2.checksum_valid(), "checksum after profile repack");
    // K single-bin chunks → K profile points read back.
    assert_eq!(rf2.profile(1).expect("profile back").point_count(), k);
    // The MS2 neighbor is intact.
    let base = RawFile::open(sample("small2.RAW")).unwrap();
    assert_eq!(base.centroid_peaks(2).len(), rf2.centroid_peaks(2).len());
    assert!(rf2.scan_filter(1).is_some());
}

/// Batch repack of many scans in one rebuild must equal doing them one-by-one, and
/// read back correctly.
#[test]
fn repack_many_grows_multiple_scans() {
    use thermorawfile::ScanEdit;
    let base = RawFile::open(sample("small2.RAW")).unwrap();
    // Pick the first four MS2/centroid scans (skip the interleaved MS1 profile scans).
    let targets: Vec<u32> = (base.first_scan..=base.last_scan)
        .filter(|&s| !base.centroid_peaks(s).is_empty())
        .take(4)
        .collect();
    assert_eq!(targets.len(), 4, "need 4 centroid scans");
    let untouched = (base.first_scan..=base.last_scan)
        .find(|&s| !base.centroid_peaks(s).is_empty() && !targets.contains(&s))
        .unwrap();
    // Distinct over-budget payloads per scan.
    let payloads: Vec<Vec<(f64, f32)>> = targets
        .iter()
        .map(|&s| {
            (0..600 + s as usize * 7)
                .map(|i| (200.0 + i as f64 * 0.3, 50.0 + i as f32))
                .collect()
        })
        .collect();
    let edits: Vec<ScanEdit> = targets
        .iter()
        .zip(&payloads)
        .map(|(&scan, p)| ScanEdit::Centroids { scan, peaks: p })
        .collect();

    // Batch: one rebuild.
    let mut a = RawFile::open(sample("small2.RAW")).unwrap();
    a.repack_many(&edits).unwrap();

    // Sequential: N splices, same edits in scan order.
    let mut b = RawFile::open(sample("small2.RAW")).unwrap();
    for (&s, p) in targets.iter().zip(&payloads) {
        b.repack_centroids(s, p).unwrap();
    }

    // Equivalence: the batch rebuild is byte-identical to the sequential result.
    assert_eq!(a.bytes.len(), b.bytes.len(), "batch vs sequential size differs");
    assert!(a.bytes == b.bytes, "batch rebuild != sequential repack (bytes differ)");

    // Read back: grown scans have their new counts; an untouched scan is intact.
    let out = format!("{}/repack_many.raw", env!("CARGO_TARGET_TMPDIR"));
    a.save(&out).unwrap();
    let rf = RawFile::open(&out).unwrap();
    assert!(rf.checksum_valid());
    for (k, &s) in targets.iter().enumerate() {
        assert_eq!(rf.centroid_peaks(s).len(), payloads[k].len(), "scan {s} count");
    }
    assert_eq!(
        base.centroid_peaks(untouched).len(),
        rf.centroid_peaks(untouched).len(),
        "untouched scan {untouched}"
    );
}

/// Batch repack must handle an UNORDERED edit list mixing grow and shrink.
#[test]
fn repack_many_unordered_with_shrink() {
    use thermorawfile::ScanEdit;
    let base = RawFile::open(sample("small2.RAW")).unwrap();
    let c: Vec<u32> = (base.first_scan..=base.last_scan)
        .filter(|&s| !base.centroid_peaks(s).is_empty())
        .take(4)
        .collect();
    // Edits NOT in scan order; mix grows and shrinks.
    let payloads: Vec<(u32, Vec<(f64, f32)>)> = vec![
        (c[3], (0..2000).map(|i| (200.0 + i as f64 * 0.1, 10.0 + i as f32)).collect()),
        (c[1], (0..5).map(|i| (300.0 + i as f64, 10.0 + i as f32)).collect()),
        (c[0], (0..1500).map(|i| (250.0 + i as f64 * 0.1, 5.0 + i as f32)).collect()),
        (c[2], (0..3).map(|i| (400.0 + i as f64, 9.0 + i as f32)).collect()),
    ];
    let edits: Vec<ScanEdit> = payloads
        .iter()
        .map(|(s, p)| ScanEdit::Centroids { scan: *s, peaks: p })
        .collect();

    let mut rf = RawFile::open(sample("small2.RAW")).unwrap();
    rf.repack_many(&edits).unwrap();
    let out = format!("{}/repack_unord.raw", env!("CARGO_TARGET_TMPDIR"));
    rf.save(&out).unwrap();
    let r = RawFile::open(&out).unwrap();
    assert!(r.checksum_valid());
    for (s, p) in &payloads {
        assert_eq!(r.centroid_peaks(*s).len(), p.len(), "scan {s} count");
    }
    // The remaining scans are untouched.
    for s in base.first_scan..=base.last_scan {
        if payloads.iter().any(|(t, _)| t == &s) {
            continue;
        }
        assert_eq!(base.centroid_peaks(s).len(), r.centroid_peaks(s).len(), "untouched {s}");
    }
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
