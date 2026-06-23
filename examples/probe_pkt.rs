//! Dump the first MS2 scan's centroid-packet header to settle wide vs narrow.
use thermorawfile::RawFile;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let rf = RawFile::open(&a[1]).expect("open");
    let mut ms2=0u32; let mut idx=0usize;
    for s in 1..=rf.scan_count() as u32 {
        if let Some(e)=rf.scan_event(s){ if e.ms_order>=2 { ms2=s; idx=(s-rf.first_scan) as usize; break; }}
    }
    let pkt=(rf.data_addr + rf.index[idx].offset) as usize;
    let u32at=|o:usize| u32::from_le_bytes(rf.bytes[o..o+4].try_into().unwrap());
    let profile_size=u32at(pkt+4);
    let peaklist_words=u32at(pkt+8);
    let count=u32at(pkt+40);
    println!("MS2 scan {ms2}: profile_size={profile_size} peaklist_words={peaklist_words} count={count}");
    if count>0 {
        let n=count as u64; let w=peaklist_words as u64;
        println!("  narrow eqn (words==1+2*count): {}  -> {}", 1+2*n, w==1+2*n);
        println!("  wide   eqn (words==1+3*count): {}  -> {}", 1+3*n, w==1+3*n);
        println!("  VERDICT: {}", if w==1+3*n {"WIDE (f64 m/z) — current narrow detect is a BUG"} else if w==1+2*n {"NARROW (f32 m/z) — current behavior correct"} else {"NEITHER — malformed/unknown layout"});
    }
}
