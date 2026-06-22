//! `probe_kernel <in.raw> <scan>` — inject ONE fragment into the largest empty profile gap of
//! `scan` via overlay_profile, and dump the frequency-grid bins before/after. Ground-truth check
//! of the deposited peak shape (delta spike vs Gaussian), with no reader/diff in the loop.
//!
//! Calibration note: this uses the FIRST scan's FTMS calibration (`scantrailer_addr + 4`), because
//! locating an arbitrary scan's event offset needs the variable-length scan-event walk. The
//! deposited *shape* — the thing this probe validates — is calibration-independent (deposit and
//! readout share the same calib). Only the reported absolute m/z assumes `<scan>` shares the first
//! scan's calibration; use the first scan for an exact m/z readout.
use thermorawfile::RawFile;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let scan: u32 = a[2].parse().unwrap();
    let mut rf = RawFile::open(&a[1]).unwrap();
    let cal = rf
        .calibration_at_event(rf.scantrailer_addr as usize + 4)
        .unwrap();
    let prof = rf.profile(scan).expect("scan has no FTMS profile");

    // largest empty gap between chunks -> a clean injection site
    let mut best = (0u32, prof.chunks[0].first_bin + prof.chunks[0].signal.len() as u32 + 50);
    for w in prof.chunks.windows(2) {
        let end = w[0].first_bin + w[0].signal.len() as u32;
        let start = w[1].first_bin;
        if start > end && start - end > best.0 {
            best = (start - end, (end + start) / 2);
        }
    }
    let mid_bin = best.1;
    let target = prof.mz_of_bin(mid_bin, &cal);
    println!("scan {scan}: injecting one fragment at m/z {target:.4} (bin {mid_bin}, gap {} bins)", best.0);

    let dump = |rf: &RawFile, label: &str| {
        let p = rf.profile(scan).unwrap();
        let mut apex = 0f32;
        for ch in &p.chunks {
            for &v in &ch.signal {
                if v > apex {
                    apex = v;
                }
            }
        }
        println!("\n{label}:");
        for ch in &p.chunks {
            for (j, &v) in ch.signal.iter().enumerate() {
                let bin = ch.first_bin + j as u32;
                if (bin as i64 - mid_bin as i64).abs() <= 8 {
                    let bar = "#".repeat((40.0 * v as f64 / apex.max(1.0) as f64) as usize);
                    println!("   bin {bin}  mz {:.4}  {:14.1} {bar}", p.mz_of_bin(bin, &cal), v);
                }
            }
        }
    };
    dump(&rf, "BEFORE (should be empty in this window)");
    rf.overlay_profile(scan, &[(target, 1.0e6f32)], &cal).unwrap();
    dump(&rf, "AFTER overlay_profile (the deposited peak shape)");
}
