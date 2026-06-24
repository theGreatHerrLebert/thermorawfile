use std::time::Instant;
use thermorawfile::{RawFile, ScanEdit};
fn main() {
    let src = "/home/administrator/thermo-raw-spike/data/astral_singlecell_L3.raw";
    let out = std::env::args().nth(1).expect("out");
    let mut rf = RawFile::open(src).expect("open");
    // Grow 200 MS2 centroid scans to ~2x their peaks in ONE rebuild.
    let targets: Vec<u32> = (rf.first_scan..=rf.last_scan)
        .filter(|&s| !rf.centroid_peaks(s).is_empty())
        .take(200).collect();
    let payloads: Vec<Vec<(f64,f32)>> = targets.iter().map(|&s| {
        let n = rf.centroid_peaks(s).len();
        (0..n*2).map(|i| (200.0 + i as f64*0.05, 100.0 + i as f32)).collect()
    }).collect();
    let edits: Vec<ScanEdit> = targets.iter().zip(&payloads)
        .map(|(&scan,p)| ScanEdit::Centroids{scan,peaks:p}).collect();
    let t = Instant::now();
    rf.repack_many(&edits).expect("repack_many");
    let dt = t.elapsed();
    rf.save(&out).expect("save");
    println!("repack_many grew {} scans in ONE rebuild in {:?}", targets.len(), dt);
    println!("check scan {} -> {} peaks (was {})", targets[0],
             RawFile::open(&out).unwrap().centroid_peaks(targets[0]).len(), payloads[0].len());
}
