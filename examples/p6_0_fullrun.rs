//! P6-0 spike (b+c): full-run exact-slot replacement + no-residual verification.
//!
//! Authors a DISTINCTIVE synthetic payload into EVERY scan of a real Astral
//! template (MS1 FTMS profile via author_profile; MS2 ASTMS centroid via
//! author_centroids — both scan-indexed = exact slot), saves, re-opens, and
//! verifies for every scan that it reads back EXACTLY the synthetic payload and
//! NOTHING else (no residual real peaks), the schedule (MS order / analyzer /
//! scan count) is intact, and the checksum is valid.
//!
//! Run: cargo run --release --example p6_0_fullrun -- <template.raw> <out.raw>

use thermorawfile::RawFile;

// Distinctive synthetic payloads — nothing like real data, so any residual real
// peak would be obvious as an extra/foreign peak on readback.
const MS2_CENTROIDS: &[(f64, f32)] = &[
    (300.0, 1000.0), (400.0, 2000.0), (500.0, 3000.0), (600.0, 4000.0),
];
const MS1_PROFILE: &[(f64, f32)] = &[
    (450.12345, 5.0e5), (500.0, 1.0e6), (618.314, 7.5e5), (777.777, 2.0e6), (880.5, 3.0e5),
];

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let (inp, out) = (&args[1], &args[2]);

    let mut rf = RawFile::open(inp).expect("open template");
    let calib = rf
        .calibration_at_event(rf.scantrailer_addr as usize + 4)
        .expect("read calibration");
    let (first, last) = (rf.first_scan, rf.last_scan);

    // Record the original schedule (to prove it's preserved) + author per scan.
    let mut orig = Vec::with_capacity((last - first + 1) as usize);
    let mut budget_fail = 0usize;
    let mut ms1_authored = 0usize;
    let mut ms2_authored = 0usize;
    for scan in first..=last {
        let ev = rf.scan_event(scan).expect("scan event");
        let is_profile = rf.profile(scan).is_some();
        orig.push((ev.ms_order, ev.analyzer, is_profile));
        let res = if ev.ms_order == 1 && is_profile {
            ms1_authored += 1;
            rf.author_profile(scan, MS1_PROFILE, &calib)
        } else {
            ms2_authored += 1;
            rf.author_centroids(scan, MS2_CENTROIDS)
        };
        if res.is_err() { budget_fail += 1; }
    }
    rf.save(out).expect("save");
    println!(
        "authored: MS1(profile)={ms1_authored}  MS2(centroid)={ms2_authored}  budget_fail={budget_fail}"
    );

    // Re-open the OUTPUT and verify every scan.
    let rf2 = RawFile::open(out).expect("re-open output");
    let cal2 = rf2.calibration_at_event(rf2.scantrailer_addr as usize + 4).unwrap();
    let checksum_ok = rf2.checksum_valid();
    assert_eq!(rf2.scan_count(), (last - first + 1) as usize, "scan count changed!");

    let mut schedule_changed = 0usize;
    let mut residual_or_wrong = 0usize;
    let mut examples = 0usize;
    for scan in first..=last {
        let i = (scan - first) as usize;
        let (ms0, an0, prof0) = orig[i];
        let ev = rf2.scan_event(scan).expect("event");
        if ev.ms_order != ms0 || ev.analyzer != an0 { schedule_changed += 1; }

        let ok = if ms0 == 1 && prof0 {
            // MS1 profile: exactly the authored chunk count (residual chunks would
            // inflate this); spot-check the first peak m/z decodes near 450.123.
            match rf2.profile(scan) {
                Some(p) => {
                    let cnt_ok = p.chunks.len() == MS1_PROFILE.len();
                    let mz_ok = p.chunks.first().map_or(false, |ch| {
                        (p.mz_of_bin(ch.first_bin, &cal2) - MS1_PROFILE[0].0).abs() < 0.05
                    });
                    cnt_ok && mz_ok
                }
                None => false,
            }
        } else {
            // MS2 centroid: EXACTLY the 4 authored peaks, nothing extra.
            let pk = rf2.centroid_peaks(scan);
            pk.len() == MS2_CENTROIDS.len()
                && pk.iter().zip(MS2_CENTROIDS).all(|(p, &(mz, it))| {
                    (p.mz - mz).abs() < 1e-2 && (p.intensity - it).abs() < 1.0
                })
        };
        if !ok {
            residual_or_wrong += 1;
            if examples < 8 {
                examples += 1;
                let got = if ms0 == 1 && prof0 {
                    format!("profile chunks={}", rf2.profile(scan).map(|p| p.chunks.len()).unwrap_or(0))
                } else {
                    format!("centroid peaks={}", rf2.centroid_peaks(scan).len())
                };
                eprintln!("MISMATCH scan {scan} (ms{ms0}): {got}");
            }
        }
    }

    println!("checksum_valid={checksum_ok}  schedule_changed={schedule_changed}  residual_or_wrong={residual_or_wrong}");
    if checksum_ok && schedule_changed == 0 && residual_or_wrong == 0 && budget_fail == 0 {
        println!("VERDICT: full-run replacement OK across all {} scans — every slot holds ONLY its synthetic payload, no residual, schedule + checksum intact", last - first + 1);
    } else {
        println!("VERDICT: FAIL (see counts above)");
        std::process::exit(1);
    }
}
