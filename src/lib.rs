//! Pure-Rust reader for Thermo Finnigan `.raw` files (no .NET / RawFileReader DLL).
//!
//! Binary layout ported from `unthermo` (Apache-2.0, Pieter Kelchtermans /
//! proteinspector) and corrected/extended against real rev-66 files using the
//! Thermo RawFileReader as an oracle. Notable correction vs. unthermo: rev-66
//! FTMS centroid peaks are `{ f64 m/z, f32 intensity }` (12 bytes), not
//! `{ f32, f32 }` — *except* the Astral analyzer (ASTMS), whose centroid peaks
//! are `{ f32 m/z, f32 intensity }` (8 bytes). The width is selected per scan
//! from the peak-list word count (see [`peak_is_wide`]).
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

/// A decoded FTMS profile: a frequency grid plus contiguous signal chunks.
///
/// The profile is stored as a sparse set of `chunks` over a uniform frequency
/// grid (`first_value` + bin·`step`). Converting a bin to m/z needs the
/// per-scan frequency→m/z calibration (not yet ported), so this struct exposes
/// the grid verbatim — enough to rewrite intensities in place
/// ([`RawFile::set_profile_intensities`]) on a template's real m/z grid.
#[derive(Clone, Debug)]
pub struct Profile {
    pub first_value: f64,
    pub step: f64,
    /// Total number of bins in the (sparse) grid.
    pub nbins: u32,
    pub chunks: Vec<ProfileChunk>,
}

impl Profile {
    /// Total number of stored signal points across all chunks.
    pub fn point_count(&self) -> usize {
        self.chunks.iter().map(|c| c.signal.len()).sum()
    }
}

/// One contiguous run of profile signal points starting at `first_bin`.
#[derive(Clone, Debug)]
pub struct ProfileChunk {
    pub first_bin: u32,
    pub fudge: f32,
    pub signal: Vec<f32>,
}

struct MsRunHeader {
    first_scan: u32,
    last_scan: u32,
    scan_index_addr: u64,
    data_addr: u64,
    scantrailer_addr: u64,
    scanparams_addr: u64,
}

// Scan-event field offsets within a fixed-size event record (rev 66), decoded
// empirically against RawFileReader (unthermo's variable-length v66 layout is
// wrong). Events are a contiguous fixed-stride array in [scantrailer+4, scanparams).
const EV_MS_ORDER: usize = 6; // preamble: 1 = MS1, 2 = MS2
const EV_ANALYZER: usize = 40; // 0 = ITMS, 4 = FTMS
const EV_ISO_CENTER: usize = 140; // f64 precursor / isolation-window center m/z
const EV_ISO_WIDTH: usize = 148; // f64 isolation width
const EV_COLLISION_ENERGY: usize = 156; // f64 collision energy

/// The acquisition descriptor for one scan: MS order, analyzer, and (for MS2)
/// the quadrupole isolation window + collision energy.
#[derive(Clone, Copy, Debug)]
pub struct ScanEvent {
    pub ms_order: u8,
    pub analyzer: u8,
    pub isolation_center: f64,
    pub isolation_width: f64,
    pub collision_energy: f64,
}

/// A parsed Thermo `.raw` file held entirely in memory.
pub struct RawFile {
    pub bytes: Vec<u8>,
    pub version: u32,
    pub first_scan: u32,
    pub last_scan: u32,
    pub scan_index_addr: u64,
    pub data_addr: u64,
    pub scantrailer_addr: u64,
    pub scanparams_addr: u64,
    /// Fixed stride of a scan-event record (bytes); 0 if it could not be derived.
    pub scan_event_size: usize,
    pub index: Vec<ScanIndexEntry>,
}

fn err(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// Width of a centroid peak record, selected per scan.
///
/// The packet's `peaklist_size` is a count of 4-byte words including the leading
/// `u32` peak count, so `(peaklist_size - 1) / count` is the words per peak:
/// FTMS-style packets (e.g. Orbitrap Velos) use 3 words = 12 bytes
/// `{ f64 m/z, f32 intensity }`; the Astral analyzer (ASTMS) uses 2 words =
/// 8 bytes `{ f32 m/z, f32 intensity }`. Returns `true` for the wide (12-byte)
/// form. Defaults to wide when the count is unknown/inconsistent.
fn peak_is_wide(peaklist_size: u32, count: u32) -> bool {
    if count == 0 {
        return true;
    }
    peaklist_size.saturating_sub(1) / count >= 3
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

        // Scan events are a fixed-stride array in [scantrailer+4, scanparams).
        let region = (ms.scanparams_addr).saturating_sub(ms.scantrailer_addr + 4) as usize;
        let scan_event_size = if n > 0 && region >= n && region % n == 0 {
            region / n
        } else {
            0
        };

        Ok(RawFile {
            bytes,
            version,
            first_scan: ms.first_scan,
            last_scan: ms.last_scan,
            scan_index_addr: ms.scan_index_addr,
            data_addr: ms.data_addr,
            scantrailer_addr: ms.scantrailer_addr,
            scanparams_addr: ms.scanparams_addr,
            scan_event_size,
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
            let wide = peak_is_wide(peaklist_size, count);
            peaks.reserve(count as usize);
            for _ in 0..count {
                let mz = if wide { c.f64() } else { c.f32() as f64 };
                let intensity = c.f32();
                peaks.push(Peak { mz, intensity });
            }
        }
        peaks
    }

    /// Decode the FTMS profile for `scan` (1-based), or `None` if the scan has
    /// no profile (centroid-only, e.g. ASTMS MS2). The profile is the chunked
    /// frequency-grid signal that precedes the centroid label list in the packet.
    pub fn profile(&self, scan: u32) -> Option<Profile> {
        if scan < self.first_scan || scan > self.last_scan {
            return None;
        }
        let e = &self.index[(scan - self.first_scan) as usize];
        let pos = (self.data_addr + e.offset) as usize;
        let mut c = Cur::new(&self.bytes, pos);
        let _unknown1 = c.u32();
        let profile_size = c.u32();
        let _peaklist_size = c.u32();
        let layout = c.u32();
        c.skip(16);
        let _low = c.f32();
        let _high = c.f32();
        if profile_size == 0 {
            return None;
        }
        let first_value = c.f64();
        let step = c.f64();
        let peak_count = c.u32();
        let nbins = c.u32();
        let mut chunks = Vec::with_capacity(peak_count as usize);
        for _ in 0..peak_count {
            let first_bin = c.u32();
            let cn = c.u32();
            // Fudge is present only when the packet layout flag is non-zero.
            let fudge = if layout > 0 { c.f32() } else { 0.0 };
            let mut signal = Vec::with_capacity(cn as usize);
            for _ in 0..cn {
                signal.push(c.f32());
            }
            chunks.push(ProfileChunk {
                first_bin,
                fudge,
                signal,
            });
        }
        Some(Profile {
            first_value,
            step,
            nbins,
            chunks,
        })
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

// ---------------------------------------------------------------------------
// Writer (template-mutation): rewrite scan data in an existing file, then fix
// the index stats and the integrity checksum. Same-count peak rewrites need no
// offset rebuild; variable counts are TODO.
// ---------------------------------------------------------------------------

fn put_u32(b: &mut [u8], off: usize, v: u32) {
    b[off..off + 4].copy_from_slice(&v.to_le_bytes());
}
fn put_f32(b: &mut [u8], off: usize, v: f32) {
    b[off..off + 4].copy_from_slice(&v.to_le_bytes());
}
fn put_f64(b: &mut [u8], off: usize, v: f64) {
    b[off..off + 8].copy_from_slice(&v.to_le_bytes());
}

impl RawFile {
    /// Overwrite a scan's centroid peak list in place. The new list must have the
    /// **same number of peaks** as the current one (variable counts require a scan-
    /// index offset rebuild — TODO). Recomputes the scan-index stats (TIC, base
    /// peak, m/z range) and the packet-header m/z range so the file stays
    /// internally consistent. Call [`RawFile::save`] afterwards to fix the checksum.
    pub fn set_centroid_peaks(&mut self, scan: u32, peaks: &[Peak]) -> io::Result<()> {
        if scan < self.first_scan || scan > self.last_scan {
            return Err(err("scan out of range"));
        }
        let idx = (scan - self.first_scan) as usize;
        let entry = self.index[idx].clone();
        let pkt = (self.data_addr + entry.offset) as usize;

        let profile_size = u32::from_le_bytes(self.bytes[pkt + 4..pkt + 8].try_into().unwrap());
        let peaklist_size = u32::from_le_bytes(self.bytes[pkt + 8..pkt + 12].try_into().unwrap());
        if profile_size != 0 || peaklist_size == 0 {
            return Err(err("not a centroid-only scan (profile rewrite is TODO)"));
        }
        let peaklist_size = u32::from_le_bytes(self.bytes[pkt + 8..pkt + 12].try_into().unwrap());
        let count = u32::from_le_bytes(self.bytes[pkt + 40..pkt + 44].try_into().unwrap()) as usize;
        if peaks.len() != count {
            return Err(err(
                "peak count must equal the existing count (variable-count rewrite is TODO)",
            ));
        }

        // Peak record width is analyzer-dependent: 12 bytes { f64 m/z, f32 int }
        // for FTMS, 8 bytes { f32 m/z, f32 int } for the Astral analyzer (ASTMS).
        let wide = peak_is_wide(peaklist_size, count as u32);
        let stride = if wide { 12 } else { 8 };
        let peaks_off = pkt + 44;
        let mut tic = 0f64;
        let mut base_mz = 0f64;
        let mut base_int = 0f64;
        let mut low_mz = f64::INFINITY;
        let mut high_mz = f64::NEG_INFINITY;
        for (i, p) in peaks.iter().enumerate() {
            // Stats must reflect the value actually stored, so for the narrow
            // form use the f32-rounded m/z (what the reader will read back).
            let mz_stored = if wide { p.mz } else { (p.mz as f32) as f64 };
            if wide {
                put_f64(&mut self.bytes, peaks_off + i * stride, p.mz);
                put_f32(&mut self.bytes, peaks_off + i * stride + 8, p.intensity);
            } else {
                put_f32(&mut self.bytes, peaks_off + i * stride, p.mz as f32);
                put_f32(&mut self.bytes, peaks_off + i * stride + 4, p.intensity);
            }
            tic += p.intensity as f64;
            if (p.intensity as f64) > base_int {
                base_int = p.intensity as f64;
                base_mz = mz_stored;
            }
            low_mz = low_mz.min(mz_stored);
            high_mz = high_mz.max(mz_stored);
        }
        if peaks.is_empty() {
            low_mz = 0.0;
            high_mz = 0.0;
        }

        // Packet header m/z range (f32 @ +32 / +36).
        put_f32(&mut self.bytes, pkt + 32, low_mz as f32);
        put_f32(&mut self.bytes, pkt + 36, high_mz as f32);

        // Scan-index entry stats (f64): TIC @ +32, base-int @ +40, base-mz @ +48,
        // low-mz @ +56, high-mz @ +64.
        let ea = self.scan_index_addr as usize + idx * scan_index_entry_size(self.version);
        put_f64(&mut self.bytes, ea + 32, tic);
        put_f64(&mut self.bytes, ea + 40, base_int);
        put_f64(&mut self.bytes, ea + 48, base_mz);
        put_f64(&mut self.bytes, ea + 56, low_mz);
        put_f64(&mut self.bytes, ea + 64, high_mz);

        // Refresh the cached entry.
        self.index[idx] = ScanIndexEntry {
            total_current: tic,
            base_mz,
            low_mz,
            high_mz,
            ..entry
        };
        Ok(())
    }

    /// Overwrite an FTMS profile's signal intensities in place, in chunk order.
    ///
    /// `signal` must have the same total length as the existing profile
    /// ([`Profile::point_count`]); the frequency grid (chunk first-bins / step /
    /// m/z calibration) is preserved, so this rewrites intensities onto the
    /// template's real m/z axis — the basis for emitting a simulated MS1 onto a
    /// real Astral/Orbitrap grid. Recomputes the scan-index TIC and base
    /// intensity; the m/z range and base-peak m/z are left unchanged (the grid
    /// is unchanged, and base-peak m/z needs the frequency→m/z calibration,
    /// which is a follow-up). Call [`RawFile::save`] to fix the checksum.
    pub fn set_profile_intensities(&mut self, scan: u32, signal: &[f32]) -> io::Result<()> {
        if scan < self.first_scan || scan > self.last_scan {
            return Err(err("scan out of range"));
        }
        let idx = (scan - self.first_scan) as usize;
        let entry = self.index[idx].clone();
        let pkt = (self.data_addr + entry.offset) as usize;

        let profile_size = u32::from_le_bytes(self.bytes[pkt + 4..pkt + 8].try_into().unwrap());
        if profile_size == 0 {
            return Err(err("scan has no profile (centroid-only)"));
        }
        let layout = u32::from_le_bytes(self.bytes[pkt + 12..pkt + 16].try_into().unwrap());

        // Walk the chunk headers to collect the byte offset of every signal value.
        let mut c = Cur::new(&self.bytes, pkt + 40);
        let _first_value = c.f64();
        let _step = c.f64();
        let peak_count = c.u32();
        let _nbins = c.u32();
        let mut offsets: Vec<usize> = Vec::new();
        for _ in 0..peak_count {
            let _first_bin = c.u32();
            let cn = c.u32();
            if layout > 0 {
                c.skip(4); // fudge
            }
            for _ in 0..cn {
                offsets.push(c.p);
                c.skip(4);
            }
        }
        if signal.len() != offsets.len() {
            return Err(err(
                "signal length must equal the existing profile point count (Profile::point_count)",
            ));
        }

        let mut tic = 0f64;
        let mut base_int = 0f64;
        for (&off, &v) in offsets.iter().zip(signal) {
            put_f32(&mut self.bytes, off, v);
            tic += v as f64;
            base_int = base_int.max(v as f64);
        }

        // Scan-index entry: TIC @ +32, base-int @ +40. base-mz/low-mz/high-mz
        // are left as-is (grid unchanged; base-peak m/z needs the calibration).
        let ea = self.scan_index_addr as usize + idx * scan_index_entry_size(self.version);
        put_f64(&mut self.bytes, ea + 32, tic);
        put_f64(&mut self.bytes, ea + 40, base_int);
        self.index[idx] = ScanIndexEntry {
            total_current: tic,
            ..entry
        };
        Ok(())
    }

    /// Recompute and write the Adler-32 integrity checksum into the header.
    pub fn recompute_checksum(&mut self) {
        let crc = compute_checksum(&self.bytes);
        put_u32(&mut self.bytes, CHECKSUM_OFFSET, crc);
    }

    /// Fix the checksum and write the file to disk.
    pub fn save<P: AsRef<Path>>(&mut self, path: P) -> io::Result<()> {
        self.recompute_checksum();
        std::fs::write(path, &self.bytes)
    }

    fn scan_event_offset(&self, scan: u32) -> Option<usize> {
        if self.scan_event_size == 0 || scan < self.first_scan || scan > self.last_scan {
            return None;
        }
        Some(
            self.scantrailer_addr as usize
                + 4
                + (scan - self.first_scan) as usize * self.scan_event_size,
        )
    }

    /// Read the acquisition descriptor (MS order, analyzer, isolation, CE) for `scan`.
    pub fn scan_event(&self, scan: u32) -> Option<ScanEvent> {
        let o = self.scan_event_offset(scan)?;
        let f64at = |off: usize| f64::from_le_bytes(self.bytes[off..off + 8].try_into().unwrap());
        Some(ScanEvent {
            ms_order: self.bytes[o + EV_MS_ORDER],
            analyzer: self.bytes[o + EV_ANALYZER],
            isolation_center: f64at(o + EV_ISO_CENTER),
            isolation_width: f64at(o + EV_ISO_WIDTH),
            collision_energy: f64at(o + EV_COLLISION_ENERGY),
        })
    }

    /// Author an MS2 isolation window: set the precursor / window-center m/z,
    /// isolation width, and collision energy for `scan`. Call [`RawFile::save`]
    /// afterwards to fix the checksum.
    pub fn set_isolation(
        &mut self,
        scan: u32,
        center: f64,
        width: f64,
        collision_energy: f64,
    ) -> io::Result<()> {
        let o = self
            .scan_event_offset(scan)
            .ok_or_else(|| err("no scan event (unknown event stride or scan out of range)"))?;
        put_f64(&mut self.bytes, o + EV_ISO_CENTER, center);
        put_f64(&mut self.bytes, o + EV_ISO_WIDTH, width);
        put_f64(&mut self.bytes, o + EV_COLLISION_ENERGY, collision_energy);
        Ok(())
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
    let scanparams_addr = c.u64();
    MsRunHeader {
        first_scan,
        last_scan,
        scan_index_addr,
        data_addr,
        scantrailer_addr,
        scanparams_addr,
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
