//! `rawdump <file.raw> [scan]` — dump structure, or a scan's centroid peaks.

use std::process::exit;
use thermorawfile::RawFile;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: rawdump <file.raw> [scan]");
        exit(2);
    }
    let rf = match RawFile::open(&args[1]) {
        Ok(rf) => rf,
        Err(e) => {
            eprintln!("error: {e}");
            exit(1);
        }
    };

    println!("version          = {}", rf.version);
    println!("scans            = {}..{} ({} total)", rf.first_scan, rf.last_scan, rf.scan_count());
    println!("scan_index_addr  = {:#x}", rf.scan_index_addr);
    println!("data_addr        = {:#x}", rf.data_addr);
    println!(
        "checksum         = stored {:#010x}, computed {:#010x}  [{}]",
        rf.stored_checksum(),
        rf.compute_checksum(),
        if rf.checksum_valid() { "VALID" } else { "MISMATCH" }
    );

    if let Some(s) = args.get(2).and_then(|s| s.parse::<u32>().ok()) {
        let peaks = rf.centroid_peaks(s);
        println!("\nscan {s}: {} centroid peaks", peaks.len());
        for (i, p) in peaks.iter().take(10).enumerate() {
            println!("  [{i}] m/z={:.4}  int={:.1}", p.mz, p.intensity);
        }
        if peaks.len() > 10 {
            let i = peaks.len() - 1;
            println!("  ... [{i}] m/z={:.4}  int={:.1}", peaks[i].mz, peaks[i].intensity);
        }
    }
}
