use thermorawfile::RawFile;
fn main() {
    for f in [concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/small2.RAW"),
              "/home/administrator/thermo-raw-spike/data/astral_singlecell_L3.raw"] {
        match RawFile::open(f) {
            Ok(rf) => println!("{}: {} controller(s)", f, rf.controller_dir.len()),
            Err(e) => println!("{}: open err {e}", f),
        }
    }
}
