# Alignment Compression: mpox-2000 xz vs zstd

The level-1 zstd alignment preserves the decoded FASTA exactly, increases compressed size from 431,188 bytes to 24,910,268 bytes, and does not produce a measurable end-to-end single-thread speedup. The full scaling results are in [mpox-2000.md](mpox-2000.md) and [mpox-2000-zstd.md](mpox-2000-zstd.md).

## Inputs

| Format | File | Size | Decoded SHA-256 |
| --- | --- | ---: | --- |
| xz | `data/mpox/clade-ii/2000/aln.fasta.xz` | 431,188 bytes | `0be5e75cb1b2c3d01c7d0c53528d3bb16493febcc5765459805677c6eeb5640e` |
| zstd level 1 | `data/mpox/clade-ii/2000/aln.fasta.zst` | 24,910,268 bytes | `0be5e75cb1b2c3d01c7d0c53528d3bb16493febcc5765459805677c6eeb5640e` |

The zstd file was produced with `xz -dc aln.fasta.xz | zstd -1 > aln.fasta.zst`. Level 1 is zstd's fast-compression preset. The resulting file is 57.8x larger than the xz input.

## Methodology

Both reports use the same five workloads, thread counts of 1, 2, 4, 8, and 16, three measured hyperfine runs after one warmup, a second duplicate run per configuration, and separate `/usr/bin/time` peak-RSS measurements. Only the alignment supplied to `ancestral`, `optimize`, and `timetree` differs. `mugration` and `clock` do not read the alignment and serve as environmental controls.

The xz report used commit `3f1820a5ee5f443691624e00ab4aad4e4012a049`; the zstd report used `7e8a78ca6d068dafba5f31946cc25becde0217ec`. No TreeTime production source changed between these commits. The zstd report records identical SHA-256 hashes for its two executable copies.

## Wall-clock comparison

Negative deltas favor zstd. Values are means from separate benchmark sessions, so deltas include session-level scheduling and system-load variation.

| Workload | Threads | xz | zstd | Delta |
| --- | ---: | ---: | ---: | ---: |
| **Ancestral** | 1 | 5.830 s | 5.738 s | -1.6% |
|  | 2 | 5.518 s | 4.999 s | -9.4% |
|  | 4 | 5.012 s | 4.610 s | -8.0% |
|  | 8 | 4.681 s | 4.720 s | +0.8% |
|  | 16 | 4.627 s | 4.854 s | +4.9% |
| **Optimize** | 1 | 5.897 s | 5.855 s | -0.7% |
|  | 2 | 4.899 s | 4.942 s | +0.9% |
|  | 4 | 4.347 s | 3.851 s | -11.4% |
|  | 8 | 4.967 s | 4.556 s | -8.3% |
|  | 16 | 4.054 s | 4.024 s | -0.7% |
| **Timetree** | 1 | 34.947 s | 35.732 s | +2.2% |
|  | 2 | 21.171 s | 21.098 s | -0.3% |
|  | 4 | 14.079 s | 13.673 s | -2.9% |
|  | 8 | 11.011 s | 12.529 s | +13.8% |
|  | 16 | 9.748 s | 14.220 s | +45.9% |

At one thread, all three differences are within 2.2% and point in mixed directions. Higher-thread deltas are also inconsistent: ancestral and optimize improve at some counts, while timetree regresses sharply at 8 and 16 threads. The no-alignment controls vary between the sessions as well, including mugration and sub-100 ms clock timings. The data therefore support no causal end-to-end speedup from zstd in this run.

## Interpretation

[profiling-findings.md](profiling-findings.md) attributes 18.8% of one-thread ancestral self-time to FASTA reading and xz decompression. That creates an upper bound on the benefit available from changing compression alone. Serial sparse-partition construction, output, page management, and other orchestration remain, while only 25.6% of runtime is parallel Fitch computation. Per-depth-level fork/join barriers continue to constrain parallel scaling regardless of the alignment codec.

Level-1 zstd is therefore a poor trade for this dataset under the measured workflow: it consumes 57.8x more storage without a detectable end-to-end runtime improvement. The codec change does preserve all alignment-dependent outputs at both tested thread extremes.
