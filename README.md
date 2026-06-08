# thermorawfile

A **pure-Rust** reader (and, soon, writer) for **Thermo Finnigan `.raw`** files —
no .NET runtime, no Thermo `RawFileReader` DLL, no Windows. Cross-platform, zero
dependencies.

## Why

Every existing way to touch a Thermo `.raw` file goes through Thermo's closed
`RawFileReader` .NET assembly (that's what ProteoWizard, ThermoRawFileParser,
`mzdata`'s `thermo` feature, etc. all call). That makes them **read-only** and
runtime-heavy. This crate parses the format natively, and is being extended into
the **only open-source `.raw` _writer_** — so synthetic / simulated acquisitions
can be emitted in a format that Thermo's own tools accept.

## Status

- ✅ Native reader for file revision ≥ 64 (Orbitrap-era): file header / version,
  run-header walk, scan index, **centroid peak lists**.
- ✅ Integrity checksum (Adler-32) compute + verify — matches Thermo's reader.
- 🚧 Writer (in progress): peak rewrite proven; scan-index rebuild + index-stat
  recompute next.
- 🚧 TODO: FTMS profile packets, file revisions < 64, more instrument variants.

## Usage

```rust
use thermorawfile::RawFile;

let rf = RawFile::open("run.raw")?;
println!("rev {} — {} scans", rf.version, rf.scan_count());
assert!(rf.checksum_valid());
for p in rf.centroid_peaks(2) {
    println!("{:.4}  {:.1}", p.mz, p.intensity);
}
```

CLI:

```
cargo run --release --bin rawdump -- run.raw 2
```

## Format notes (the bits that were not obvious)

- The container is the proprietary **Finnigan** format (magic `0xA101` +
  `"Finnigan"` UTF-16), **not** OLE2/MCDF. Structural layout follows
  [`unthermo`](https://bitbucket.org/proteinspector/ms) (Apache-2.0).
- **Integrity checksum:** a 4-byte little-endian **Adler-32 at file offset 148**,
  computed over `file[0 : min(len, 10 MiB)]` with that field zeroed, seeded with
  `0` (i.e. `zlib.adler32(buf, 0)`). Thermo's reader rejects any file whose bytes
  don't match — so a writer must recompute it.
- **Centroid peaks (rev 66):** 12 bytes each — `float64 m/z` followed by
  `float32 intensity`. (Older revisions use `float32` m/z; `unthermo` assumes the
  8-byte form, which is why rev-66 peaks decode wrong there.)

## Tests

`cargo test` validates the reader against real sample files in `tests/data/`
(from the Apache-2.0 [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser)
test corpus), cross-checked against values from Thermo's own `RawFileReader`.

## License

Dual-licensed under **MIT OR Apache-2.0**. Binary layout knowledge derives from
the Apache-2.0 `unthermo` project; see [`NOTICE`](NOTICE).
