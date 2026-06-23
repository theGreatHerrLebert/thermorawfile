//! Self-describing "GenericRecord" streams in Thermo `.raw` (trailer scan-parameters and
//! status / instrument logs). A [`GenericDataHeader`] declares a schema of typed, labelled
//! fields; each following record is a fixed-stride row decoded against that schema. This is
//! the substrate for per-scan trailer metadata (AGC, ion-injection time, FT resolution,
//! charge state, FAIMS, NCE, …).
//!
//! Ported from **OpenTFRaw** (`generic_data.rs` / `bytes.rs` / `reader.rs`), Apache-2.0,
//! Copyright Sigilweaver Holdings LLC — itself clean-room reverse-engineered from public
//! PRIDE deposits (see NOTICE). Adapted to thermorawfile's in-memory byte model and
//! `io::Result` error type; parsing behaviour is preserved. No Thermo SDK / proprietary code.

use std::io::{self, ErrorKind};

fn eof(offset: usize, needed: usize) -> io::Error {
    io::Error::new(ErrorKind::UnexpectedEof, format!("need {needed} bytes at offset {offset}"))
}
fn bad_utf16(offset: usize) -> io::Error {
    io::Error::new(ErrorKind::InvalidData, format!("invalid UTF-16 at offset {offset}"))
}

/// Fallible cursor over an in-memory byte slice. Every read is bounds-checked — no
/// `slice[a..b].try_into().unwrap()` chains (codex: no unchecked slicing on metadata paths).
pub(crate) struct SliceReader<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> SliceReader<'a> {
    pub(crate) fn new(buf: &'a [u8]) -> Self {
        Self { buf, pos: 0 }
    }
    pub(crate) fn position(&self) -> usize {
        self.pos
    }
    pub(crate) fn len(&self) -> usize {
        self.buf.len()
    }
    pub(crate) fn seek_to(&mut self, off: usize) -> io::Result<()> {
        if off > self.buf.len() {
            return Err(eof(off, 0));
        }
        self.pos = off;
        Ok(())
    }
    fn take(&mut self, n: usize) -> io::Result<&'a [u8]> {
        let s = self.buf.get(self.pos..self.pos + n).ok_or_else(|| eof(self.pos, n))?;
        self.pos += n;
        Ok(s)
    }
    fn array<const N: usize>(&mut self) -> io::Result<[u8; N]> {
        let p = self.pos;
        self.take(N)?.try_into().map_err(|_| eof(p, N))
    }
    pub(crate) fn read_bytes(&mut self, n: usize) -> io::Result<Vec<u8>> {
        Ok(self.take(n)?.to_vec())
    }
    fn read_u8(&mut self) -> io::Result<u8> {
        Ok(self.take(1)?[0])
    }
    fn read_i8(&mut self) -> io::Result<i8> {
        Ok(self.take(1)?[0] as i8)
    }
    fn read_u16(&mut self) -> io::Result<u16> {
        Ok(u16::from_le_bytes(self.array()?))
    }
    fn read_i16(&mut self) -> io::Result<i16> {
        Ok(i16::from_le_bytes(self.array()?))
    }
    fn read_u32(&mut self) -> io::Result<u32> {
        Ok(u32::from_le_bytes(self.array()?))
    }
    fn read_i32(&mut self) -> io::Result<i32> {
        Ok(i32::from_le_bytes(self.array()?))
    }
    fn read_f32(&mut self) -> io::Result<f32> {
        Ok(f32::from_le_bytes(self.array()?))
    }
    fn read_f64(&mut self) -> io::Result<f64> {
        Ok(f64::from_le_bytes(self.array()?))
    }

    /// Fixed-width UTF-16LE string of `byte_len` bytes, null-stripped.
    fn read_utf16_fixed(&mut self, byte_len: usize) -> io::Result<String> {
        let p = self.pos;
        if byte_len % 2 != 0 {
            return Err(bad_utf16(p));
        }
        let raw = self.take(byte_len)?;
        let units: Vec<u16> = raw.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]])).collect();
        let end = units.iter().position(|&u| u == 0).unwrap_or(units.len());
        String::from_utf16(&units[..end]).map_err(|_| bad_utf16(p))
    }

    /// PascalStringWin32: u32 character count, then that many UTF-16LE code units.
    fn read_pascal_string(&mut self) -> io::Result<String> {
        let p = self.pos;
        let char_count = self.read_u32()? as usize;
        if char_count == 0 {
            return Ok(String::new());
        }
        let byte_len = char_count.checked_mul(2).ok_or_else(|| bad_utf16(p))?;
        let raw = self.take(byte_len)?;
        let units: Vec<u16> = raw.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]])).collect();
        let end = units.iter().position(|&u| u == 0).unwrap_or(units.len());
        String::from_utf16(&units[..end]).map_err(|_| bad_utf16(p))
    }
}

/// Field value type codes used in a GenericDataHeader schema.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenericType {
    Gap,
    Int8,
    Bool,
    BoolYesNo,
    BoolOnOff,
    UInt8,
    Int16,
    UInt16,
    Int32,
    UInt32,
    Float32,
    Float64,
    AsciiString,
    WideString,
}

impl GenericType {
    fn from_u32(v: u32) -> Option<Self> {
        Some(match v {
            0x0 => Self::Gap,
            0x1 => Self::Int8,
            0x2 => Self::Bool,
            0x3 => Self::BoolYesNo,
            0x4 => Self::BoolOnOff,
            0x5 => Self::UInt8,
            0x6 => Self::Int16,
            0x7 => Self::UInt16,
            0x8 => Self::Int32,
            0x9 => Self::UInt32,
            0xA => Self::Float32,
            0xB => Self::Float64,
            0xC => Self::AsciiString,
            0xD => Self::WideString,
            _ => return None,
        })
    }
}

/// One field descriptor within a [`GenericDataHeader`].
#[derive(Debug, Clone)]
pub struct GenericDataDescriptor {
    pub field_type: GenericType,
    pub length: u32,
    pub label: String,
}

/// Self-describing schema for a GenericRecord stream.
#[derive(Debug)]
pub struct GenericDataHeader {
    pub fields: Vec<GenericDataDescriptor>,
}

impl GenericDataHeader {
    /// Try to read a header at the current position. Returns `Ok(None)` (restoring the
    /// position) if the bytes do not look like a valid schema — implausible field count,
    /// unknown type code, or non-string label. Fails closed on truncation.
    pub(crate) fn try_read(r: &mut SliceReader) -> io::Result<Option<Self>> {
        let saved = r.position();
        let n = r.read_u32()?;
        if !(2..=500).contains(&n) {
            r.seek_to(saved)?;
            return Ok(None);
        }
        let mut fields = Vec::with_capacity(n as usize);
        for _ in 0..n {
            let type_code = r.read_u32()?;
            let field_type = match GenericType::from_u32(type_code) {
                Some(t) => t,
                None => {
                    r.seek_to(saved)?;
                    return Ok(None);
                }
            };
            let length = r.read_u32()?;
            // Peek the label's character count; a sane label is short.
            let label_start = r.position();
            let char_count = r.read_u32()?;
            if char_count > 200 {
                r.seek_to(saved)?;
                return Ok(None);
            }
            r.seek_to(label_start)?;
            let label = match r.read_pascal_string() {
                Ok(s) => s,
                Err(e) if e.kind() == ErrorKind::InvalidData => {
                    r.seek_to(saved)?;
                    return Ok(None);
                }
                Err(e) => return Err(e),
            };
            if label.len() > 200 {
                r.seek_to(saved)?;
                return Ok(None);
            }
            fields.push(GenericDataDescriptor { field_type, length, label });
        }
        let hdr = Self { fields };
        if !hdr.looks_meaningful() {
            r.seek_to(saved)?;
            return Ok(None);
        }
        Ok(Some(hdr))
    }

    fn looks_meaningful(&self) -> bool {
        let named = self.fields.iter().filter(|f| !f.label.is_empty()).count();
        named >= 2 && self.fixed_record_size() > 0
    }

    /// Fixed on-disk byte size of one record under this schema.
    pub(crate) fn fixed_record_size(&self) -> usize {
        self.fields
            .iter()
            .map(|f| match f.field_type {
                // Gap = padding of `length` bytes (OpenTFRaw treats it as 0-width, which
                // misaligns any schema with a non-zero-length gap — codex review).
                GenericType::Gap => f.length as usize,
                GenericType::Int8
                | GenericType::Bool
                | GenericType::BoolYesNo
                | GenericType::BoolOnOff
                | GenericType::UInt8 => 1,
                GenericType::Int16 | GenericType::UInt16 => 2,
                GenericType::Int32 | GenericType::UInt32 | GenericType::Float32 => 4,
                GenericType::Float64 => 8,
                GenericType::AsciiString => f.length as usize,
                GenericType::WideString => f.length as usize * 2,
            })
            .sum()
    }

    /// Locate a plausible header within a bounded forward window. The v64+ error-log region
    /// holds padding of indeterminate size before the scan-parameters schema, so the schema
    /// is found by signature scan. On success `r` is positioned at the first record.
    pub(crate) fn find_forward(
        r: &mut SliceReader,
        max_scan: usize,
        expected_record_size: Option<usize>,
    ) -> io::Result<Option<Self>> {
        let start = r.position();
        let cap = max_scan.min(4 * 1024 * 1024).min(r.len().saturating_sub(start));
        let buf = r.read_bytes(cap)?; // owned copy — avoids per-candidate borrows of `r`
        for pass in 0..2 {
            let mut off = 0usize;
            while off + 4 <= buf.len() {
                let n = u32::from_le_bytes([buf[off], buf[off + 1], buf[off + 2], buf[off + 3]]);
                if (2..=500).contains(&n) {
                    let mut sub = SliceReader::new(&buf[off..]);
                    if let Some(hdr) = Self::try_read(&mut sub)? {
                        let size_ok =
                            pass == 1 || expected_record_size.is_none_or(|w| hdr.fixed_record_size() == w);
                        if size_ok {
                            r.seek_to(start + off)?;
                            let _ = Self::try_read(r)?; // consume the header → `r` at first record
                            return Ok(Some(hdr));
                        }
                    }
                }
                off += 2;
            }
            if expected_record_size.is_none() {
                break;
            }
        }
        r.seek_to(start)?;
        Ok(None)
    }
}

/// A typed value decoded from a generic record.
#[derive(Debug, Clone)]
pub enum GenericValue {
    Gap,
    Int8(i8),
    Bool(bool),
    UInt8(u8),
    Int16(i16),
    UInt16(u16),
    Int32(i32),
    UInt32(u32),
    Float32(f32),
    Float64(f64),
    String(String),
}

impl GenericValue {
    /// Numeric value as f64, where meaningful.
    pub fn as_f64(&self) -> Option<f64> {
        Some(match self {
            Self::Float64(v) => *v,
            Self::Float32(v) => *v as f64,
            Self::Int32(v) => *v as f64,
            Self::UInt32(v) => *v as f64,
            Self::Int16(v) => *v as f64,
            Self::UInt16(v) => *v as f64,
            Self::Int8(v) => *v as f64,
            Self::UInt8(v) => *v as f64,
            _ => return None,
        })
    }
}

/// One record decoded against a [`GenericDataHeader`].
#[derive(Debug)]
pub struct GenericRecord {
    pub values: Vec<(String, GenericValue)>,
}

impl GenericRecord {
    pub(crate) fn read(r: &mut SliceReader, header: &GenericDataHeader) -> io::Result<Self> {
        let mut values = Vec::with_capacity(header.fields.len());
        for desc in &header.fields {
            let value = match desc.field_type {
                GenericType::Gap => {
                    let _ = r.read_bytes(desc.length as usize)?; // consume the padding
                    GenericValue::Gap
                }
                GenericType::Int8 => GenericValue::Int8(r.read_i8()?),
                GenericType::Bool | GenericType::BoolYesNo | GenericType::BoolOnOff => {
                    GenericValue::Bool(r.read_u8()? != 0)
                }
                GenericType::UInt8 => GenericValue::UInt8(r.read_u8()?),
                GenericType::Int16 => GenericValue::Int16(r.read_i16()?),
                GenericType::UInt16 => GenericValue::UInt16(r.read_u16()?),
                GenericType::Int32 => GenericValue::Int32(r.read_i32()?),
                GenericType::UInt32 => GenericValue::UInt32(r.read_u32()?),
                GenericType::Float32 => GenericValue::Float32(r.read_f32()?),
                GenericType::Float64 => GenericValue::Float64(r.read_f64()?),
                GenericType::AsciiString => {
                    let s = if desc.length > 0 {
                        let bytes = r.read_bytes(desc.length as usize)?;
                        let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
                        String::from_utf8_lossy(&bytes[..end]).into_owned()
                    } else {
                        String::new()
                    };
                    GenericValue::String(s)
                }
                GenericType::WideString => {
                    let s = if desc.length > 0 {
                        r.read_utf16_fixed(desc.length as usize * 2)?
                    } else {
                        String::new()
                    };
                    GenericValue::String(s)
                }
            };
            values.push((desc.label.clone(), value));
        }
        Ok(Self { values })
    }

    /// Look up a field value by exact label.
    pub fn get(&self, label: &str) -> Option<&GenericValue> {
        self.values.iter().find(|(l, _)| l == label).map(|(_, v)| v)
    }
    pub fn get_f64(&self, label: &str) -> Option<f64> {
        self.get(label)?.as_f64()
    }
    pub fn get_i32(&self, label: &str) -> Option<i32> {
        match self.get(label)? {
            GenericValue::Int32(v) => Some(*v),
            GenericValue::Int16(v) => Some(*v as i32),
            GenericValue::Int8(v) => Some(*v as i32),
            GenericValue::UInt32(v) => i32::try_from(*v).ok(),
            GenericValue::UInt16(v) => Some(*v as i32),
            GenericValue::UInt8(v) => Some(*v as i32),
            _ => None,
        }
    }
    pub fn get_bool(&self, label: &str) -> Option<bool> {
        match self.get(label)? {
            GenericValue::Bool(v) => Some(*v),
            _ => None,
        }
    }
    pub fn get_string(&self, label: &str) -> Option<&str> {
        match self.get(label)? {
            GenericValue::String(v) => Some(v.as_str()),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pascal(s: &str) -> Vec<u8> {
        let units: Vec<u16> = s.encode_utf16().collect();
        let mut out = (units.len() as u32).to_le_bytes().to_vec();
        for u in units {
            out.extend_from_slice(&u.to_le_bytes());
        }
        out
    }
    fn field(type_code: u32, length: u32, label: &str) -> Vec<u8> {
        let mut out = type_code.to_le_bytes().to_vec();
        out.extend_from_slice(&length.to_le_bytes());
        out.extend(pascal(label));
        out
    }

    #[test]
    fn parses_synthetic_schema_and_record() {
        let mut hdr = 2u32.to_le_bytes().to_vec();
        hdr.extend(field(0xB, 0, "Ion Injection Time (ms):")); // Float64
        hdr.extend(field(0x8, 0, "Charge State:")); // Int32
        let mut r = SliceReader::new(&hdr);
        let header = GenericDataHeader::try_read(&mut r).unwrap().expect("valid header");
        assert_eq!(header.fields.len(), 2);
        assert_eq!(header.fixed_record_size(), 12); // f64(8) + i32(4)

        let mut rec = 12.5f64.to_le_bytes().to_vec();
        rec.extend_from_slice(&3i32.to_le_bytes());
        let mut rr = SliceReader::new(&rec);
        let record = GenericRecord::read(&mut rr, &header).unwrap();
        assert_eq!(record.get_f64("Ion Injection Time (ms):"), Some(12.5));
        assert_eq!(record.get_i32("Charge State:"), Some(3));
        assert!(record.get_string("missing").is_none());
    }

    #[test]
    fn find_forward_locates_schema_after_padding() {
        let mut blob = vec![0xAAu8; 36]; // junk / error-log gap (even → 2-byte aligned schema)
        let here = blob.len();
        blob.extend(2u32.to_le_bytes());
        blob.extend(field(0xB, 0, "Resolution:"));
        blob.extend(field(0x8, 0, "Charge State:"));
        blob.extend(99.0f64.to_le_bytes()); // first record begins here
        let mut r = SliceReader::new(&blob);
        let header = GenericDataHeader::find_forward(&mut r, 4096, Some(12)).unwrap().expect("found");
        assert_eq!(header.fields.len(), 2);
        assert_eq!(r.position(), here + 4 + (12 + "Resolution:".len() * 2) + (12 + "Charge State:".len() * 2));
        assert_eq!(r.read_f64().unwrap(), 99.0); // positioned at the first record
    }

    #[test]
    fn rejects_bogus_header() {
        let bytes = 9999u32.to_le_bytes();
        let mut r = SliceReader::new(&bytes);
        assert!(GenericDataHeader::try_read(&mut r).unwrap().is_none());
        assert_eq!(r.position(), 0); // position restored
    }

    #[test]
    fn truncation_errors_not_panics() {
        let bytes = [0x02u8, 0x00]; // claims a field count but is truncated
        let mut r = SliceReader::new(&bytes);
        assert!(r.read_u32().is_err());
    }

    #[test]
    fn gap_with_length_is_consumed() {
        let mut hdr = 3u32.to_le_bytes().to_vec();
        hdr.extend(field(0xB, 0, "A:")); // Float64
        hdr.extend(field(0x0, 4, "g:")); // Gap, 4 bytes
        hdr.extend(field(0x8, 0, "B:")); // Int32
        let mut r = SliceReader::new(&hdr);
        let header = GenericDataHeader::try_read(&mut r).unwrap().expect("valid header");
        assert_eq!(header.fixed_record_size(), 8 + 4 + 4);

        let mut rec = 1.5f64.to_le_bytes().to_vec();
        rec.extend_from_slice(&[0xFF; 4]); // gap padding
        rec.extend_from_slice(&7i32.to_le_bytes());
        let mut rr = SliceReader::new(&rec);
        let record = GenericRecord::read(&mut rr, &header).unwrap();
        assert_eq!(record.get_f64("A:"), Some(1.5));
        assert_eq!(record.get_i32("B:"), Some(7)); // 0xFFFFFFFF garbage if the gap is not skipped
    }

    #[test]
    fn get_i32_reads_uint8() {
        let mut hdr = 2u32.to_le_bytes().to_vec();
        hdr.extend(field(0x5, 0, "Charge:")); // UInt8
        hdr.extend(field(0x5, 0, "N:")); // UInt8
        let mut r = SliceReader::new(&hdr);
        let header = GenericDataHeader::try_read(&mut r).unwrap().expect("valid header");
        let mut rr = SliceReader::new(&[4u8, 9u8]);
        let rec = GenericRecord::read(&mut rr, &header).unwrap();
        assert_eq!(rec.get_i32("Charge:"), Some(4));
        assert_eq!(rec.get_i32("N:"), Some(9));
    }
}
