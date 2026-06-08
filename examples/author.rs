//! `author <in.raw> <out.raw> <scan>` — read a file, overwrite one scan's centroid
//! peaks with a synthetic ramp, and save (recomputing index stats + checksum).
//! Used to prove that a Thermo .raw written entirely by Rust is read back correctly
//! by Thermo's own RawFileReader.

use thermorawfile::{Peak, RawFile};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 4 {
        eprintln!("usage: author <in.raw> <out.raw> <scan>");
        std::process::exit(2);
    }
    let scan: u32 = a[3].parse().expect("scan number");
    let mut rf = RawFile::open(&a[1]).expect("open input");
    let n = rf.centroid_peaks(scan).len();
    let peaks: Vec<Peak> = (0..n)
        .map(|i| Peak {
            mz: 250.0 + i as f64 * 0.5,
            intensity: 700.0 + i as f32 * 4.0,
        })
        .collect();
    rf.set_centroid_peaks(scan, &peaks).expect("rewrite peaks");
    rf.save(&a[2]).expect("save output");
    println!("wrote {n} synthetic peaks into scan {scan} -> {}", a[2]);
}
