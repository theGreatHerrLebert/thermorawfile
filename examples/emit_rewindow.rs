use thermorawfile::RawFile;
fn main() {
    let src = "/home/administrator/thermo-raw-spike/data/RD139_Narrow_Fusion_DIA.raw";
    let out = std::env::args().nth(1).unwrap();
    let mut rf = RawFile::open(src).unwrap();
    // Re-window scan 2: 354@8Th  ->  new center 400.0, width 4.0, CE 30.
    rf.set_isolation(2, 400.0, 4.0, 30.0).unwrap();
    rf.save(&out).unwrap();
    println!("re-windowed scan 2 -> center 400.0 width 4.0 CE 30; wrote {out}");
}
