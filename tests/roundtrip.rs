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

/// Non-gated: scan_event_ranges decodes on the fixed-stride fixture, and
/// rewindow_in_place skips MS1 + honors a `None` callback (re-windows only the one scan).
#[test]
fn scan_event_ranges_and_rewindow_on_fixture() {
    let mut rf = RawFile::open(sample("small2.RAW")).unwrap();
    let ms2 = (rf.first_scan..=rf.last_scan)
        .find(|&s| rf.scan_event(s).is_some_and(|e| e.ms_order >= 2))
        .expect("a MS2 scan in the fixture");
    let ranges = rf.scan_event_ranges(ms2).expect("range block decodes on fixed-stride");
    assert!(
        !ranges.is_empty() && ranges.iter().all(|&(lo, hi)| lo > 0.0 && hi > lo),
        "plausible fragment scan range, got {ranges:?}"
    );

    let mut closure_saw_ms1 = false;
    let n = rf
        .rewindow_in_place(|s, e| {
            if e.ms_order == 1 {
                closure_saw_ms1 = true; // must never fire — the loop filters MS1 out
            }
            (s == ms2).then_some((500.0, 3.0, 20.0)) // re-window only this scan; None elsewhere
        })
        .unwrap();
    assert!(!closure_saw_ms1, "rewindow_in_place must skip MS1 (callback never sees one)");
    assert_eq!(n, 1, "only the one scan re-windowed; None left the rest");
    let ev = rf.scan_event(ms2).unwrap();
    assert!((ev.isolation_center - 500.0).abs() < 1e-6 && (ev.isolation_width - 3.0).abs() < 1e-6);
}

/// Tier-2 3a: re-window all MS2 scans in place (here: halve every isolation width) and
/// read the new windows back. Gate on a DIA template: `TIMSIM_VARLEN_DIA_TEMPLATE=<dia .raw>`.
#[test]
fn rewindow_in_place_roundtrip() {
    let p = match std::env::var("TIMSIM_VARLEN_DIA_TEMPLATE") {
        Ok(p) => p,
        Err(_) => {
            eprintln!("SKIP rewindow_in_place_roundtrip: set TIMSIM_VARLEN_DIA_TEMPLATE=<dia .raw>");
            return;
        }
    };
    let before = RawFile::open(&p).unwrap();
    // Record the first MS2 scan's original window + a neighbor's.
    let first_ms2 = (before.first_scan..=before.last_scan)
        .find(|&s| before.scan_event(s).is_some_and(|e| e.ms_order >= 2))
        .expect("a MS2 scan");
    let orig = before.scan_event(first_ms2).unwrap();

    let mut rf = RawFile::open(&p).unwrap();
    // Halve each MS2 window's width; keep center + CE.
    let n = rf
        .rewindow_in_place(|_s, e| Some((e.isolation_center, e.isolation_width / 2.0, e.collision_energy)))
        .unwrap();
    assert!(n > 0, "expected MS2 scans to re-window");

    let out = format!("{}/rewindowed.raw", env!("CARGO_TARGET_TMPDIR"));
    rf.save(&out).unwrap();
    let r = RawFile::open(&out).unwrap();
    assert!(r.checksum_valid(), "checksum after re-window");
    let now = r.scan_event(first_ms2).unwrap();
    assert!((now.isolation_center - orig.isolation_center).abs() < 1e-6, "center preserved");
    assert!((now.isolation_width - orig.isolation_width / 2.0).abs() < 1e-6, "width halved");
    assert_eq!(r.scan_event_ranges(first_ms2), before.scan_event_ranges(first_ms2), "scan range unchanged");
    eprintln!("rewindow OK: {n} MS2 re-windowed; width {} -> {}", orig.isolation_width, now.isolation_width);
}

/// Variable-length scan events (Orbitrap Fusion-class, no fixed stride) must decode to
/// the correct per-scan ms-order + isolation windows. Gate on such a file:
/// `TIMSIM_VARLEN_TEMPLATE=<fusion .raw>`.
#[test]
fn variable_length_scan_events_decode() {
    let p = match std::env::var("TIMSIM_VARLEN_TEMPLATE") {
        Ok(p) => p,
        Err(_) => {
            eprintln!("SKIP variable_length_scan_events_decode: set TIMSIM_VARLEN_TEMPLATE=<fusion-class .raw>");
            return;
        }
    };
    let rf = RawFile::open(&p).unwrap();
    assert!(rf.has_scan_events(), "variable-length scan events must decode to an offset table");
    let mut ms1 = 0usize;
    let mut ms2_with_iso = 0usize;
    for s in rf.first_scan..=rf.last_scan {
        if let Some(e) = rf.scan_event(s) {
            if e.ms_order == 1 {
                ms1 += 1;
            } else if e.ms_order >= 2 && e.isolation_center > 0.0 && e.isolation_width > 0.0 {
                ms2_with_iso += 1;
            }
        }
    }
    // The old fixed-stride-only code mislabeled these as all-MS2 with no isolation.
    assert!(ms1 > 0, "expected decoded MS1 scans, got none");
    assert!(ms2_with_iso > 0, "expected MS2 scans with decoded isolation windows, got none");
    eprintln!("variable-length OK: {ms1} MS1, {ms2_with_iso} MS2 with isolation windows");
}

/// A controller-directory RunHeaderAddr pointing past EOF (malformed or foreign-layout
/// file) must produce a graceful Err — NOT a slice panic, which would abort the whole
/// process when called across the PyO3 boundary.
#[test]
fn open_errs_not_panics_on_out_of_range_runheader_addr() {
    let mut bytes = std::fs::read(sample("small2.RAW")).unwrap();
    assert!(RawFile::from_bytes(bytes.clone()).is_ok(), "fixture opens clean as-is");
    // The single controller's RunHeaderAddr (u64) for this fixture sits at offset 2432.
    let off = 2432usize;
    assert_eq!(
        u64::from_le_bytes(bytes[off..off + 8].try_into().unwrap()),
        2_071_234,
        "fixture RunHeaderAddr offset moved; update the test"
    );
    bytes[off..off + 8].copy_from_slice(&u64::MAX.to_le_bytes());
    assert!(
        RawFile::from_bytes(bytes).is_err(),
        "out-of-range RunHeaderAddr must error gracefully, not panic"
    );
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

/// DIRECT byte-preservation proof (single-controller fixture). The gap region
/// (data_end..scan_index) holds the MS run header + instrument method + tune + logs,
/// and the tail (scan-trailer..EOF) holds scan events + scan params. A repack
/// RELOCATES both and bumps ONLY the run-header's 64-bit section pointers — the
/// method/tune/log/event/param PAYLOAD must be byte-identical. We prove it by asserting
/// the only bytes that differ in the gap fall inside the run-header section-pointer
/// block, and the tail is fully identical. Independent of any reader's tolerance.
#[test]
fn repack_preserves_gap_section_payload() {
    fn data_end(rf: &RawFile) -> usize {
        rf.data_addr as usize
            + (0..rf.index.len()).map(|j| rf.index[j].data_packet_size as usize).sum::<usize>()
    }
    let before = RawFile::open(sample("small2.RAW")).unwrap();
    assert_eq!(before.controller_dir.len(), 1, "this byte-proof assumes one controller");
    let de0 = data_end(&before);
    // Run-header offset within the gap, and its relocatable 64-bit pointer block
    // (scan_index @ +7408 .. scanparams @ +7456, each u64). Bytes here SHOULD change.
    let rh_off = before.ms_runheader_addr as usize - de0;
    let ptr_block = (rh_off + 7408)..(rh_off + 7456 + 8);
    let tail_orig = before.bytes[before.scantrailer_addr as usize..].to_vec();

    let mut rf = RawFile::open(sample("small2.RAW")).unwrap();
    let grown: Vec<(f64, f32)> = (0..600).map(|i| (250.0 + i as f64 * 0.5, 700.0 + i as f32)).collect();
    rf.repack_centroids(2, &grown).unwrap();

    let g0 = &before.bytes[de0..before.scan_index_addr as usize];
    let g1 = &rf.bytes[data_end(&rf)..rf.scan_index_addr as usize];
    assert_eq!(g0.len(), g1.len(), "gap length changed (should only relocate)");
    let mut changed_outside = 0;
    for i in 0..g0.len() {
        if g0[i] != g1[i] && !ptr_block.contains(&i) {
            changed_outside += 1;
        }
    }
    assert_eq!(
        changed_outside, 0,
        "gap payload (method/tune/log) bytes changed OUTSIDE the run-header pointer block"
    );

    // The scan-trailer/params tail is relocated but never rewritten → fully identical.
    let tail_new = &rf.bytes[rf.scantrailer_addr as usize..];
    assert!(tail_orig == tail_new, "scan-trailer/params BYTES changed after repack");
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
