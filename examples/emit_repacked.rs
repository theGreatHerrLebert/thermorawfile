//! Emit a repacked .raw for RawFileReader validation.
//! Run: cargo run --release --example emit_repacked -- <in.raw> <out.raw> <scan> <mult>
use thermorawfile::RawFile;

fn main() {
    let mut a = std::env::args().skip(1);
    let src = a.next().expect("in.raw");
    let out = a.next().expect("out.raw");
    let scan: u32 = a.next().map(|s| s.parse().unwrap()).unwrap_or(2);
    let mult: usize = a.next().map(|s| s.parse().unwrap()).unwrap_or(3);

    let mut rf = RawFile::open(&src).expect("open");
    let orig = rf.centroid_peaks(scan).len();
    assert!(orig > 0, "scan {scan} has no centroid peaks (profile scan?)");

    let grown: Vec<(f64, f32)> = (0..orig * mult)
        .map(|i| (250.0 + i as f64 * 0.5, 700.0 + i as f32 * 4.0))
        .collect();
    rf.repack_centroids(scan, &grown).expect("repack");
    rf.save(&out).expect("save");
    println!("wrote {out}: scan {scan} grown {orig} -> {} peaks", grown.len());
}
