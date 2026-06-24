//! probe_format <in.raw> — report, per scan-type, the analyzer + MS2 packet format
//! (centroid vs profile) so we can scope authoring support for a new instrument.
use thermorawfile::RawFile;
use std::collections::BTreeMap;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let rf = RawFile::open(&a[1]).expect("open");
    println!("scans={} checksum_valid={}", rf.scan_count(), rf.checksum_valid());
    // analyzer code: 0=ITMS, 4=FTMS, (Astral=ASTMS, other code)
    let mut seen: BTreeMap<(u8,u8), (usize, usize, usize)> = BTreeMap::new(); // (ms_order,analyzer)->(count,with_centroids,with_profile)
    for scan in 1..=rf.scan_count() as u32 {
        let ev = match rf.scan_event(scan) { Some(e)=>e, None=>continue };
        let nc = rf.centroid_peaks(scan).len();
        let prof = rf.profile(scan).is_some();
        let e = seen.entry((ev.ms_order, ev.analyzer)).or_insert((0,0,0));
        e.0 += 1; if nc>0 { e.1+=1; } if prof { e.2+=1; }
    }
    println!("(ms_order, analyzer[0=ITMS,4=FTMS,other=ASTMS]) -> count, with_centroids, with_profile");
    for ((mso,an),(c,cn,pr)) in &seen {
        println!("  MS{} analyzer={} : {} scans | {} have centroids | {} have profile", mso, an, c, cn, pr);
    }
    // sample one MS2 scan's first few centroid m/z to confirm readability
    for scan in 1..=rf.scan_count() as u32 {
        if let Some(ev)=rf.scan_event(scan) { if ev.ms_order>=2 {
            let p=rf.centroid_peaks(scan);
            println!("sample MS2 scan {}: iso={:.2} ce={:.1} npeaks={} first_mz={:?}",
                scan, ev.isolation_center, ev.collision_energy, p.len(),
                p.iter().take(3).map(|x| (x.mz*1e4).round()/1e4).collect::<Vec<_>>());
            break;
        }}
    }
}
