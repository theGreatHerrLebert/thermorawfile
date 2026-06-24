# Repack gap-section preservation — validation record

**Question (Codex review #6):** the repack relocates the data section and shifts the
trailing "gap" sections — MS run header, instrument method, tune data, status/error
logs — plus the scan index, scan-trailer (events) and scan-params. It bumps each
*section-start* pointer but does **not** rewrite the bytes *inside* those sections. Do
any of them carry internal byte-offset pointers that would go stale when relocated?

## Two independent checks

### 1. Direct byte-range proof (authoritative) — `repack_preserves_gap_section_payload`

A pure-Rust test (no reader involved). After growing a scan on `small2.RAW`
(single-controller, rev 66):

- The gap region (`data_end .. scan_index`) is byte-identical between original and
  repacked **except** for the bytes inside the MS run-header's 64-bit section-pointer
  block (`+7408 .. +7464`) — which the repack is *supposed* to bump. Measured: of
  217,936 gap bytes, exactly **11 differ**, all inside that pointer block. The
  instrument-method / tune / log **payload is byte-identical**.
- The tail (`scan_trailer .. EOF` — scan events + scan params) is **fully
  byte-identical**: relocated, never rewritten.

This is a true byte-preservation proof for those regions: the only changes are the
section pointers we deliberately relocate; nothing inside the payloads moves or mutates.

### 2. RawFileReader cross-check (official reader, semantic projection)

A `.NET` probe `digest` mode (official Thermo `RawFileReader`) emits a diffable digest —
instrument-method streams, run header, tune count, and per-scan filter + MS order +
trailer-extra — run on original vs repacked and diffed:

- **Velos Pro** (grow 1 scan 196→588): digest **100% identical**, 0 differing lines
  (4 method streams, run header, tune, all 95 scans' trailer-extra).
- **Astral** (1 GB, **5 controllers**, batch grow 200 scans): gap sections
  (run header / method / tune) **identical**; **0 differing lines** across all
  **18,909** scans' trailer-extra + filter + MS order.

This exercises the official reader's traversal of the relocated sections end to end.

## Claim, precisely scoped (per Codex)

> For the tested files, a repack preserves the relocated gap-section / scan-event /
> scan-param **payload bytes** verbatim (proof 1), and the official RawFileReader reads
> all of its metadata projections identically (proof 2). The only mutated bytes outside
> the data section + scan index are the run-header section pointers we relocate.

We do **not** claim universal correctness for every consumer or every file shape.

## Residual gaps (not yet covered)

- The byte-range proof is on a **single-controller** fixture; the multi-controller
  (Astral) case is covered at the **reader** level (proof 2) but not yet byte-level for
  each non-MS controller's own run-header pointer blocks.
- **Cross-reader**: ProteoWizard `msconvert`/`msaccess`, Sage, DIA-NN not yet run on a
  repacked file (a stricter consumer could in principle depend on layout we haven't
  checked).
- The digest hashes RawFileReader's *projection* (FNV over `label=value`), not raw
  bytes — which is why proof 1 (direct bytes) is the authoritative one.
- A **negative control** (deliberately stale a pointer, confirm the checks fail) would
  further harden the test harness.

## Tooling

- `tests/roundtrip.rs::repack_preserves_gap_section_payload` (in-repo, runs in CI).
- `../probe` (`.NET`, not version-controlled): `dotnet probe.dll <file> digest`.
