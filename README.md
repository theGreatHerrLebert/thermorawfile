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
- **Centroid peaks (rev 66):** the record width is **per scan**, read from the
  packet header — `1+2·n` words → 8-byte (`float32 m/z`, `float32 int`),
  `1+3·n` words → 12-byte (`float64 m/z`, `float32 int`), `n` = peak count.
  `unthermo` assumes the 8-byte form, which is why rev-66 wide scans decode wrong
  there. See `centroid_record_width` and `examples/probe_pkt.rs`.

## Tests

`cargo test` validates the reader against real sample files in `tests/data/`
(from the Apache-2.0 [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser)
test corpus): the Adler-32 checksum recomputes to the stored value, and parsed
peaks are self-consistent with each scan's own structure (record width fixed by
the packet word count; `float64` m/z in range). `examples/probe_pkt.rs` reproduces
the width determination on any `.raw`.

## Provenance & reverse engineering

This crate contains **no Thermo SDK, no `RawFileReader` DLL, and no Thermo
proprietary code**. The `.raw` binary layout is reconstructed entirely from
clean-room, publicly available sources:

- [`unthermo`](https://bitbucket.org/proteinspector/ms) (Apache-2.0) — structural
  layout and field names;
- [`OpenTFRaw`](https://github.com/Sigilweaver/OpenTFRaw) (Apache-2.0) and
  `unfinnigan` — independent reverse-engineering of the v66 scan-event layout,
  preamble offsets, and frequency↔m/z calibration, both derived from public
  PRIDE deposits;
- the files' own **self-describing structure** — e.g. the centroid record width
  is fixed *per scan* by the peak-list word count, read from the packet header,
  not assumed.

The format was **validated against genuine public PRIDE deposits** (e.g. PXD060431
Orbitrap, PXD061065 Astral): the word count determines the record width per scan
and the `float64` m/z decode yields monotonic fragment m/z within each scan's own
bounds. No proprietary software is needed to obtain or verify any of this; the
determination is reproducible via `examples/probe_pkt.rs`.

## License

Dual-licensed under **MIT OR Apache-2.0**. Binary layout knowledge derives from
the Apache-2.0 `unthermo` project; see [`NOTICE`](NOTICE).
