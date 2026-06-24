//! Dump the section map of a .raw so we know what sits after the data region
//! (the repack blast radius). Run: cargo run --example layout_probe -- <file>
use thermorawfile::RawFile;

fn main() {
    let path = std::env::args().nth(1).unwrap_or_else(|| {
        format!("{}/tests/data/small2.RAW", env!("CARGO_MANIFEST_DIR"))
    });
    let rf = RawFile::open(&path).expect("open");
    let flen = std::fs::metadata(&path).unwrap().len();

    // Where does the contiguous data region end? data_addr + max(offset + size).
    let mut data_end = rf.data_addr;
    let mut min_off = u64::MAX;
    for e in &rf.index {
        let abs = rf.data_addr + e.offset;
        data_end = data_end.max(abs + e.data_packet_size as u64);
        min_off = min_off.min(e.offset);
    }

    println!("file: {}", path);
    println!("len:               {}", flen);
    println!("version:           {}", rf.version);
    println!("first/last scan:   {} .. {}", rf.first_scan, rf.last_scan);
    println!("scans:             {}", rf.index.len());
    println!("scan_event_size:   {}", rf.scan_event_size);
    println!("--- section addresses (sorted) ---");
    let mut secs: Vec<(&str, u64)> = vec![
        ("scan_index_addr", rf.scan_index_addr),
        ("data_addr",       rf.data_addr),
        ("scantrailer_addr", rf.scantrailer_addr),
        ("scanparams_addr",  rf.scanparams_addr),
    ];
    secs.sort_by_key(|s| s.1);
    for (n, a) in &secs {
        println!("  {:>16} = {:>12}   ({:+} from data_addr)", n, a, *a as i64 - rf.data_addr as i64);
    }
    println!("--- data region ---");
    println!("  data_addr        = {}", rf.data_addr);
    println!("  first entry off  = {} (min off = {})", rf.index[0].offset, min_off);
    println!("  data_end (calc)  = {}", data_end);
    println!("  bytes after data_end to EOF = {}", flen as i64 - data_end as i64);
    println!("--- ordering verdict ---");
    let after: Vec<&str> = secs.iter().filter(|(_, a)| *a > rf.data_addr).map(|(n, _)| *n).collect();
    if after.is_empty() {
        println!("  DATA IS LAST: no tracked section starts after data_addr.");
        println!("  (still verify nothing untracked — method/log/tune — sits past data_end)");
    } else {
        println!("  sections AFTER data_addr: {:?}", after);
        println!("  -> growing a packet shifts these; their pointers must move.");
    }
}
