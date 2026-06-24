**Findings**
**Medium: test claims “contiguous” but only proves monotonic.**  
In `repack_preserves_offset32_zero_and_offsets_contiguous`, the final assertion checks only:
```rust
rf.index[j].offset > rf.index[j - 1].offset
```
That does not prove packets are contiguous. A gap would still pass. To prove the claim, assert something like previous offset + previous packet size equals current offset, using the scan-index `DataPacketSize` / in-memory packet size model. The current `Offset(+72) == rf.index[j].offset` assertion is also fairly weak because the repack code writes disk from the same in-memory model; it proves serialization consistency, not independent correctness.
**Low/Medium: `cur32 != 0` is a reasonable preserve-don’t-fabricate heuristic, but it cannot distinguish “populated mirror with legitimate zero” from “structural zero”.**  
At `ea + 0`, the conditional update only happens when the existing field is non-zero:
```rust
let cur32 = ...
if cur32 != 0 {
    let off32 = u32::try_from(off)?;
    put_u32(&mut self.bytes, ea, off32);
}
```
For rev >= 64 files where `Offset32` is structurally zero, this is sound and avoids fabricating stale/unsupported mirrors. For a quirk file that populates `Offset32`, the first scan may legitimately have offset `0`, so that single entry would not be recognized as populated. Usually that is harmless if the first packet remains at offset 0, but the convention detection is per-entry, not per-file. If supporting quirk files matters, a stronger rule would be: detect whether any `Offset32` entry is non-zero before rewriting, then rewrite all entries under that mode.
The `>4GB` path is now intentionally conditional. That is correct for structural-zero rev >= 64 files: a >4GB 64-bit offset should not fail just because a legacy mirror field exists in the format. For a populated-mirror file, it still fails loudly before truncation. Shrink is also fine: non-zero mirrors get rewritten downward; structural zero remains zero.
**Low: run-header Addr32 removal looks consistent with the stated evidence, but the code now fully commits to that format assumption.**  
Removing the mirror relocation near `RH_SCAN_INDEX` et al. is correct if those preceding u32 slots are counts/flags or reserved-zero for all supported rev >= 64 files. The validated Velos Pro and Astral evidence is strong, and RawFileReader compatibility is a good practical check. The remaining risk is a rev >= 64 variant with non-zero 32-bit address mirrors consumed by some reader. This diff would leave those stale. If you want a guardrail, add a debug/test assertion over fixtures that the alleged Addr32 slots are zero or known non-address fields, rather than only documenting it.
**Correctness regression vs prior version**
The main behavioral change is intentional: zero `Offset32` fields are no longer populated, and run-header 32-bit slots are no longer relocated. That fixes the prior truncation/fabrication issue for rev >= 64 structural-zero files.
I do not see a new truncation regression in the conditional `Offset32` path. The one thing still weak is coverage: the added test proves “zeros stay zero” for one fixture and that disk `Offset(+72)` matches the in-memory index, but it does not prove true contiguity, does not exercise a populated `Offset32` fixture, and does not exercise the `>4GB` error branch.
