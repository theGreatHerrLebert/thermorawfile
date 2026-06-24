//! P6-0 spike (a): characterize a real Astral template's event manifest + packet
//! budget — the schedule the simulator must match and the per-packet capacity it
//! must fit. Read-only.
//!
//! Run: cargo run --release --example p6_0_manifest -- <template.raw>

use std::collections::BTreeMap;
use thermorawfile::RawFile;

fn analyzer_name(a: u8) -> &'static str {
    match a {
        0 => "ITMS",
        1 => "TQMS",
        2 => "SQMS",
        3 => "TOFMS",
        4 => "FTMS",
        5 => "Sector",
        7 => "ASTMS",
        _ => "?",
    }
}

fn main() {
    let path = std::env::args().nth(1).expect("usage: p6_0_manifest <template.raw>");
    let rf = RawFile::open(&path).expect("open template");
    let n = rf.scan_count();
    println!("template: {path}");
    println!("scans: {} (first={}, last={}), checksum_valid={}", n, rf.first_scan, rf.last_scan, rf.checksum_valid());

    // Per-(ms_order, analyzer, packet-type) aggregates + packet-size distribution.
    let mut by_kind: BTreeMap<(u8, u8, &str), (usize, Vec<u32>, Vec<usize>)> = BTreeMap::new();
    let mut ms1 = 0usize;
    let mut ms2 = 0usize;
    let mut rt_monotonic = true;
    let mut last_t = f64::MIN;
    let mut ce_min = f64::MAX;
    let mut ce_max = f64::MIN;
    let mut isow_min = f64::MAX;
    let mut isow_max = f64::MIN;
    let mut order_str = String::new(); // first ~60 scans' MS-order pattern (cycle)

    for scan in rf.first_scan..=rf.last_scan {
        let ev = match rf.scan_event(scan) { Some(e) => e, None => continue };
        let idx = &rf.index[(scan - rf.first_scan) as usize];
        if idx.time + 1e-9 < last_t { rt_monotonic = false; }
        last_t = idx.time;

        let is_profile = rf.profile(scan).is_some();
        let ptype = if is_profile { "profile" } else { "centroid" };
        let peaks = if is_profile {
            rf.profile(scan).map(|p| p.point_count()).unwrap_or(0)
        } else {
            rf.centroid_peaks(scan).len()
        };

        let e = by_kind.entry((ev.ms_order, ev.analyzer, ptype)).or_insert((0, Vec::new(), Vec::new()));
        e.0 += 1;
        e.1.push(idx.data_packet_size);
        e.2.push(peaks);

        if ev.ms_order == 1 { ms1 += 1; } else if ev.ms_order == 2 {
            ms2 += 1;
            ce_min = ce_min.min(ev.collision_energy);
            ce_max = ce_max.max(ev.collision_energy);
            isow_min = isow_min.min(ev.isolation_width);
            isow_max = isow_max.max(ev.isolation_width);
        }
        if (scan - rf.first_scan) < 60 { order_str.push_str(&ev.ms_order.to_string()); }
    }

    println!("MS1={ms1}  MS2={ms2}  RT_monotonic={rt_monotonic}");
    println!("MS2 isolation_width: {isow_min:.3}..{isow_max:.3}  collision_energy: {ce_min:.3}..{ce_max:.3}");
    println!("cycle pattern (first 60 scans' MS order): {order_str}");
    println!("\nper (ms_order, analyzer, packet) — count | data_packet_size bytes [min/med/max] | peaks [min/med/max]:");
    for ((ms, an, pt), (count, mut sizes, mut peaks)) in by_kind {
        sizes.sort_unstable();
        peaks.sort_unstable();
        let med = |v: &Vec<_>| if v.is_empty() { 0 } else { v[v.len()/2] } ;
        let (smin, smed, smax) = (sizes.first().copied().unwrap_or(0), med(&sizes), sizes.last().copied().unwrap_or(0));
        let (pmin, pmed, pmax) = (peaks.first().copied().unwrap_or(0), peaks[peaks.len()/2.min(peaks.len().saturating_sub(1))], peaks.last().copied().unwrap_or(0));
        println!(
            "  MS{ms} {:5} {pt:8}  n={count:6}  bytes[{smin}/{smed}/{smax}]  peaks[{pmin}/{pmed}/{pmax}]",
            analyzer_name(an)
        );
    }
}
