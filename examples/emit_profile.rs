//! Emit a profile-repacked .raw for RawFileReader validation:
//! grow the MS1 profile scan 1 of small2.RAW to ~4000 on-grid peaks.
use thermorawfile::RawFile;

fn main() {
    let out = std::env::args().nth(1).expect("usage: emit_profile <out.raw>");
    let src = format!("{}/tests/data/small2.RAW", env!("CARGO_MANIFEST_DIR"));
    let mut rf = RawFile::open(&src).expect("open");
    let cal = rf
        .calibration_at_event(rf.scantrailer_addr as usize + 4)
        .expect("MS1 calibration");
    let prof = rf.profile(1).expect("profile");
    let (fv, step, nbins) = (prof.first_value, prof.step, prof.nbins);

    let mut peaks: Vec<(f64, f32)> = Vec::new();
    let mut bin = 1u32;
    while peaks.len() < 4000 && bin < nbins {
        let mz = cal.mz(fv + bin as f64 * step);
        if let Some(f) = cal.freq(mz) {
            let rb = ((f - fv) / step).round();
            if rb >= 0.0 && rb < nbins as f64 && rb as u32 == bin {
                let n = peaks.len();
                peaks.push((mz, 1000.0 + n as f32));
            }
        }
        bin += 1;
    }
    let kk = peaks.len();
    rf.repack_profile(1, &peaks, &cal).expect("repack_profile");
    rf.save(&out).expect("save");
    println!("wrote {out}: MS1 scan 1 profile grown to {kk} peaks (was 3032 points)");
}
