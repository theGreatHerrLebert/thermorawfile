// Tier-2 step 1: map where the isolation window is (re)encoded, per MS2 scan.
use thermorawfile::RawFile;
fn main() {
    let path = std::env::args().nth(1).unwrap();
    let rf = RawFile::open(&path).unwrap();
    println!("== {} ==", path.rsplit('/').next().unwrap());
    let ms2: Vec<u32> = (rf.first_scan..=rf.last_scan)
        .filter(|&s| rf.scan_event(s).map(|e| e.ms_order >= 2).unwrap_or(false))
        .take(3).collect();
    for s in ms2 {
        let e = rf.scan_event(s).unwrap();
        let win = (e.isolation_center - e.isolation_width/2.0, e.isolation_center + e.isolation_width/2.0);
        let i = (s - rf.first_scan) as usize;
        let idx_lo = rf.index[i].low_mz; let idx_hi = rf.index[i].high_mz;
        let ranges = rf.scan_event_ranges(s);
        let filter = rf.scan_filter(s);
        let prov = rf.scan_params(s);
        println!("scan {s}:");
        println!("  reaction:    center={:.3} width={:.3} CE={:.2}  -> window [{:.3},{:.3}]", e.isolation_center, e.isolation_width, e.collision_energy, win.0, win.1);
        println!("  event ranges: {:?}", ranges);
        println!("  scan-index low/high m/z: [{:.3}, {:.3}]", idx_lo, idx_hi);
        println!("  filter: {:?}", filter);
        if let Some(p) = prov {
            println!("  trailer: monoiso_mz={:?} charge={:?} master_scan={:?}", p.monoisotopic_mz(), p.charge_state(), p.master_scan_number());
        }
    }
    // Q5: is the filter [lo-hi] / window text cached as a serialized string anywhere?
    if let Some(f) = rf.scan_filter(rf.first_scan + 1) {
        if let (Some(lb), Some(rb)) = (f.find('['), f.find(']')) {
            let needle = f[lb..=rb].as_bytes();
            let hits = rf.bytes.windows(needle.len()).filter(|w| *w == needle).count();
            // also UTF-16LE form
            let u16s: Vec<u8> = f[lb..=rb].encode_utf16().flat_map(|c| c.to_le_bytes()).collect();
            let hits16 = rf.bytes.windows(u16s.len()).filter(|w| *w == u16s.as_slice()).count();
            println!("  filter-text '{}' byte-hits: ascii={} utf16le={}", &f[lb..=rb], hits, hits16);
        }
    }
}
