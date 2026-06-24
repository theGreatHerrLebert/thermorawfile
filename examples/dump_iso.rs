//! dump_iso <in.raw> — print isolation center + CE for the first N MS2 scans.
//! Diagnostic: on a DDA/DIA run the isolation should VARY scan-to-scan; a constant
//! value across scans means the reader is misreading the scan-event isolation field.
use thermorawfile::RawFile;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let rf = RawFile::open(&a[1]).expect("open");
    let mut shown = 0;
    let mut isos: Vec<f64> = Vec::new();
    for scan in 1..=rf.scan_count() as u32 {
        if let Some(ev) = rf.scan_event(scan) {
            if ev.ms_order >= 2 {
                if shown < 25 {
                    println!("scan {:6} ms{} analyzer={} iso={:.3} ce={:.1}",
                        scan, ev.ms_order, ev.analyzer, ev.isolation_center, ev.collision_energy);
                }
                isos.push(ev.isolation_center);
                shown += 1;
            }
        }
    }
    let distinct: std::collections::BTreeSet<i64> = isos.iter().map(|x| (x*100.0) as i64).collect();
    println!("--> MS2 scans={} distinct_iso_centers={}", isos.len(), distinct.len());
}
