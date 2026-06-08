//! `author_mz <in.raw> <out.raw>` — author peaks at *arbitrary chosen m/z* into
//! an MS1 FTMS profile, using the decoded frequency↔m/z calibration.
//!
//! For each target m/z we invert the calibration (m/z → frequency → grid bin),
//! place a spike at the nearest existing grid bin, zero the rest, and save.
//! Thermo's RawFileReader should then report dominant peaks at (essentially)
//! the target m/z — proving the calibration lets us place signal at a chosen
//! m/z, not just rewrite existing intensities. (Placing on the *nearest existing*
//! bin avoids a chunk rebuild; arbitrary off-grid bins need variable-count
//! profile authoring, which is the remaining step.)

use thermorawfile::RawFile;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 3 {
        eprintln!("usage: author_mz <in.raw> <out.raw>");
        std::process::exit(2);
    }
    let scan: u32 = 1; // first MS1 scan; its event is at scantrailer_addr + 4
    let targets = [500.0_f64, 650.0, 800.0];

    let mut rf = RawFile::open(&a[1]).expect("open input");
    let prof = rf.profile(scan).expect("scan 1 has no FTMS profile");
    let calib = rf
        .calibration_at_event(rf.scantrailer_addr as usize + 4)
        .expect("read MS1 calibration");
    println!(
        "calibration: nparam={} a={:.4e} b={:.6e} c={:.6e}",
        calib.nparam, calib.a, calib.b, calib.c
    );

    // Sanity: inverse round-trips.
    for &t in &targets {
        let f = calib.freq(t).unwrap();
        println!("  m/z {t} -> f {f:.4} -> m/z {:.6} (Δ {:.1e})", calib.mz(f), (calib.mz(f) - t).abs());
    }

    // Flat list of (bin, position-in-signal-stream) for all stored points.
    let mut bins: Vec<u32> = Vec::with_capacity(prof.point_count());
    for ch in &prof.chunks {
        for j in 0..ch.signal.len() {
            bins.push(ch.first_bin + j as u32);
        }
    }
    let mut signal = vec![0.0f32; bins.len()];

    // Place a distinct spike at the existing bin nearest each target m/z.
    println!("\nplacing spikes:");
    for (k, &t) in targets.iter().enumerate() {
        let want_bin = prof.bin_of_mz(t, &calib).expect("target unreachable");
        // nearest stored bin to want_bin
        let (idx, &nearest) = bins
            .iter()
            .enumerate()
            .min_by_key(|(_, &b)| (b as i64 - want_bin).abs())
            .unwrap();
        let spike = 1.0e6 * (k as f32 + 1.0);
        signal[idx] = spike;
        let got_mz = prof.mz_of_bin(nearest, &calib);
        println!(
            "  target {t:.3} -> bin {want_bin} -> nearest stored bin {nearest} = m/z {got_mz:.4}  (spike {spike:.0})"
        );
    }

    rf.set_profile_intensities(scan, &signal).expect("write profile");
    rf.save(&a[2]).expect("save");
    println!("\nsaved -> {}  (open with the probe to confirm peaks at ~{:?})", a[2], targets);
}
