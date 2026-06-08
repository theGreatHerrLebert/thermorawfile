//! Pure-Rust reader for Thermo Finnigan `.raw` files (no .NET / RawFileReader DLL).
//!
//! Binary layout ported from `unthermo` (Apache-2.0, Pieter Kelchtermans /
//! proteinspector) and corrected/extended against real rev-66 files using the
//! Thermo RawFileReader as an oracle. Notable correction vs. unthermo: rev-66
//! centroid peaks are `{ f64 m/z, f32 intensity }` (12 bytes), not `{ f32, f32 }`.
//!
//! Scope of this foundation: structural chain + centroid peak lists + the
//! Adler-32 integrity checksum, for file revisions >= 64 (Orbitrap-era). Profile
//! (FTMS) packets and rev < 64 run-header layout are TODO.

use std::io;
use std::path::Path;

/// Size of the fixed file header that precedes the sequencer row.
pub const FILE_HEADER_SIZE: usize = 1356;
/// Offset of the 4-byte little-endian Adler-32 checksum inside the file header.
pub const CHECKSUM_OFFSET: usize = 148;
/// The checksum covers at most the first 10 MiB of the file.
pub const CHECKSUM_LIMIT: usize = 10_485_760;

/// A little-endian sequential cursor over a byte slice.
struct Cur<'a> {
    b: &'a [u8],
    p: usize,
}

impl<'a> Cur<'a> {
    fn new(b: &'a [u8], p: usize) -> Self {
        Cur { b, p }
    }
    fn u16(&mut self) -> u16 {
        let v = u16::from_le_bytes(self.b[self.p..self.p + 2].try_into().unwrap());
        self.p += 2;
        v
    }
    fn u32(&mut self) -> u32 {
        let v = u32::from_le_bytes(self.b[self.p..self.p + 4].try_into().unwrap());
        self.p += 4;
        v
    }
    fn i32(&mut self) -> i32 {
        let v = i32::from_le_bytes(self.b[self.p..self.p + 4].try_into().unwrap());
        self.p += 4;
        v
    }
    fn u64(&mut self) -> u64 {
        let v = u64::from_le_bytes(self.b[self.p..self.p + 8].try_into().unwrap());
        self.p += 8;
        v
    }
    fn f32(&mut self) -> f32 {
        let v = f32::from_le_bytes(self.b[self.p..self.p + 4].try_into().unwrap());
        self.p += 4;
        v
    }
    fn f64(&mut self) -> f64 {
        let v = f64::from_le_bytes(self.b[self.p..self.p + 8].try_into().unwrap());
        self.p += 8;
        v
    }
    fn skip(&mut self, n: usize) {
        self.p += n;
    }
    /// A Thermo PascalString: i32 length (in UTF-16 code units) + len*2 bytes.
    fn skip_pascal(&mut self) {
        let n = self.i32();
        if n > 0 {
            self.p += (n as usize) * 2;
        }
    }
}

/// One entry of the scan index (rev >= 64 layout, 88 bytes for rev 66).
#[derive(Clone, Debug)]
pub struct ScanIndexEntry {
    pub data_packet_size: u32,
    /// Offset of the scan's data packet, relative to the run header's `data_addr`.
    pub offset: u64,
    pub time: f64,
    pub total_current: f64,
    pub base_mz: f64,
    pub low_mz: f64,
    pub high_mz: f64,
}

/// A single centroid peak.
#[derive(Clone, Copy, Debug)]
pub struct Peak {
    pub mz: f64,
    pub intensity: f32,
}

struct MsRunHeader {
    first_scan: u32,
    last_scan: u32,
    scan_index_addr: u64,
    data_addr: u64,
    scantrailer_addr: u64,
}

/// A parsed Thermo `.raw` file held entirely in memory.
pub struct RawFile {
    pub bytes: Vec<u8>,
    pub version: u32,
    pub first_scan: u32,
    pub last_scan: u32,
    pub scan_index_addr: u64,
    pub data_addr: u64,
    pub index: Vec<ScanIndexEntry>,
}

fn err(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

impl RawFile {
    /// Read and parse a `.raw` file (structural chain + scan index).
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let bytes = std::fs::read(path)?;
        Self::from_bytes(bytes)
    }

    pub fn from_bytes(bytes: Vec<u8>) -> io::Result<Self> {
        if bytes.len() < FILE_HEADER_SIZE {
            return Err(err("file shorter than header"));
        }
        let version = read_version(&bytes);
        if version < 64 {
            return Err(err(
                "this foundation supports file revision >= 64 (Orbitrap-era); rev < 64 is TODO",
            ));
        }

        // Walk FileHeader -> SequencerRow -> AutoSamplerInfo -> RawFileInfo to
        // recover the run-header addresses.
        let mut c = Cur::new(&bytes, FILE_HEADER_SIZE);
        // SequencerRow (v >= 60)
        c.skip(64); // InjectionData fixed part
        for _ in 0..13 {
            c.skip_pascal();
        }
        c.skip_pascal();
        c.skip_pascal();
        c.skip_pascal();
        c.u32();
        for _ in 0..15 {
            c.skip_pascal();
        }
        // AutoSamplerInfo
        c.skip(24);
        c.skip_pascal();
        // RawFileInfo preamble
        c.u32(); // method-file-present
        c.skip(16); // 8x u16 date
        c.u32(); // unknown1
        c.u32(); // data_addr32
        let nctrl = c.u32();
        c.u32();
        c.u32();
        c.u32();
        c.skip(764); // Padding1
        c.u64(); // data_addr (64-bit)
        c.u64(); // unknown6
        let mut runheaders = Vec::with_capacity(nctrl as usize);
        for _ in 0..nctrl {
            runheaders.push(c.u64()); // RunHeaderAddr
            c.u64(); // unknown7
        }

        // The MS device is the run header whose scan-trailer address is non-zero.
        let ms = runheaders
            .iter()
            .map(|&a| read_runheader(&bytes, a as usize))
            .find(|rh| rh.scantrailer_addr != 0)
            .ok_or_else(|| err("no MS run header found"))?;

        let n = (ms.last_scan - ms.first_scan + 1) as usize;
        let entry_size = scan_index_entry_size(version);
        let mut index = Vec::with_capacity(n);
        let mut ic = Cur::new(&bytes, ms.scan_index_addr as usize);
        for _ in 0..n {
            let start = ic.p;
            ic.skip(20); // Offset32, Index, Scanevent, Scansegment, Next, Unknown1
            let data_packet_size = ic.u32();
            let time = ic.f64();
            let total_current = ic.f64();
            let _base_intensity = ic.f64();
            let base_mz = ic.f64();
            let low_mz = ic.f64();
            let high_mz = ic.f64();
            let offset = ic.u64();
            ic.p = start + entry_size;
            index.push(ScanIndexEntry {
                data_packet_size,
                offset,
                time,
                total_current,
                base_mz,
                low_mz,
                high_mz,
            });
        }

        Ok(RawFile {
            bytes,
            version,
            first_scan: ms.first_scan,
            last_scan: ms.last_scan,
            scan_index_addr: ms.scan_index_addr,
            data_addr: ms.data_addr,
            index,
        })
    }

    pub fn scan_count(&self) -> usize {
        self.index.len()
    }

    /// Read the centroid peak list for `scan` (1-based). Returns an empty vec for
    /// profile-only scans (FTMS profile decoding is TODO).
    pub fn centroid_peaks(&self, scan: u32) -> Vec<Peak> {
        if scan < self.first_scan || scan > self.last_scan {
            return Vec::new();
        }
        let e = &self.index[(scan - self.first_scan) as usize];
        let pos = (self.data_addr + e.offset) as usize;
        let mut c = Cur::new(&self.bytes, pos);
        let _unknown1 = c.u32();
        let profile_size = c.u32();
        let peaklist_size = c.u32();
        let _layout = c.u32();
        c.skip(16); // descriptor/unknown/triplet stream sizes + unknown2
        let _low = c.f32();
        let _high = c.f32();
        let mut peaks = Vec::new();
        if profile_size == 0 && peaklist_size > 0 {
            let count = c.u32();
            peaks.reserve(count as usize);
            for _ in 0..count {
                let mz = c.f64();
                let intensity = c.f32();
                peaks.push(Peak { mz, intensity });
            }
        }
        peaks
    }

    /// The checksum stored in the file header.
    pub fn stored_checksum(&self) -> u32 {
        stored_checksum(&self.bytes)
    }

    /// Recompute the checksum over the current bytes.
    pub fn compute_checksum(&self) -> u32 {
        compute_checksum(&self.bytes)
    }

    /// Whether the stored checksum matches the content (what RawFileReader verifies).
    pub fn checksum_valid(&self) -> bool {
        self.stored_checksum() == self.compute_checksum()
    }
}

fn read_version(b: &[u8]) -> u32 {
    u32::from_le_bytes(b[36..40].try_into().unwrap())
}

/// Run-header parse for rev >= 64 (64-bit addresses).
fn read_runheader(b: &[u8], addr: usize) -> MsRunHeader {
    // SampleInfo: FirstScanNumber @ +8, LastScanNumber @ +12.
    let first_scan = u32::from_le_bytes(b[addr + 8..addr + 12].try_into().unwrap());
    let last_scan = u32::from_le_bytes(b[addr + 12..addr + 16].try_into().unwrap());
    let mut c = Cur::new(b, addr);
    c.skip(592); // SampleInfo
    c.skip(6 * 520); // Filename1..6
    c.skip(16); // Unknown1, Unknown2 (f64)
    c.skip(7 * 520); // Filename7..13
    c.skip(40); // ScantrailerAddr32 .. Unknown8 (10x u32)
    let scan_index_addr = c.u64();
    let data_addr = c.u64();
    c.u64(); // InstlogAddr
    c.u64(); // ErrorlogAddr
    c.u64(); // Unknown9
    let scantrailer_addr = c.u64();
    MsRunHeader {
        first_scan,
        last_scan,
        scan_index_addr,
        data_addr,
        scantrailer_addr,
    }
}

fn scan_index_entry_size(version: u32) -> usize {
    match version {
        v if v < 64 => 72,
        64 => 80,
        _ => 88,
    }
}

/// The 4-byte little-endian Adler-32 stored at [CHECKSUM_OFFSET].
pub fn stored_checksum(bytes: &[u8]) -> u32 {
    u32::from_le_bytes(bytes[CHECKSUM_OFFSET..CHECKSUM_OFFSET + 4].try_into().unwrap())
}

/// Compute the Thermo integrity checksum: Adler-32 (seed 0) over the first
/// `min(len, 10 MiB)` bytes, with the 4-byte checksum field zeroed.
pub fn compute_checksum(bytes: &[u8]) -> u32 {
    let n = bytes.len().min(CHECKSUM_LIMIT);
    let mut buf = bytes[..n].to_vec();
    for i in CHECKSUM_OFFSET..CHECKSUM_OFFSET + 4 {
        buf[i] = 0;
    }
    adler32_seed0(&buf)
}

/// Adler-32 with a zero seed (Thermo's non-standard initialisation; matches
/// `zlib.adler32(data, 0)`).
fn adler32_seed0(data: &[u8]) -> u32 {
    const BASE: u32 = 65521;
    let mut a: u32 = 0;
    let mut b: u32 = 0;
    for &x in data {
        a = (a + x as u32) % BASE;
        b = (b + a) % BASE;
    }
    (b << 16) | a
}
