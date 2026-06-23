//! dump_scan_params <in.raw>  — per-scan trailer parameters for the first few MS2 scans.
//! Prints the full label set once (to confirm accessor labels) then the typed accessors.
use thermorawfile::RawFile;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let rf = RawFile::open(&a[1]).expect("open");
    println!("rev {}  scans {}  first_scan {}", rf.version, rf.scan_count(), rf.first_scan);
    let mut shown = 0u32;
    for s in rf.first_scan..=rf.last_scan {
        let ev = match rf.scan_event(s) {
            Some(e) => e,
            None => continue,
        };
        if ev.ms_order < 2 {
            continue;
        }
        let p = match rf.scan_params(s) {
            Some(p) => p,
            None => continue,
        };
        if shown == 0 {
            println!("--- all trailer labels (scan {s}) ---");
            for (label, val) in &p.record().values {
                println!("    {:<34} {:?}", label, val);
            }
            println!("--- typed accessors ---");
        }
        println!(
            "scan {s} ms{}: iit={:?}ms charge={:?} monoMz={:?} AGC={:?} isoW={:?}",
            ev.ms_order,
            p.ion_injection_time_ms(),
            p.charge_state(),
            p.monoisotopic_mz(),
            p.agc_target(),
            p.isolation_width_mz()
        );
        shown += 1;
        if shown >= 5 {
            break;
        }
    }
    if shown == 0 {
        println!("(no scan params found)");
    }
}
