# Dense marginal reconstruction keeps four profile arrays per node, and first-touch page faults take a fifth of its CPU time

Dense marginal reconstruction stores an L-by-s array of probabilities for every node and three more for every edge, where L is the alignment length and s the number of states. Every array is a new allocation, and the first write to each of its pages is a page fault in which the kernel zeroes the page. On a large dense run, the kernel takes about 21 percent of all CPU time, mostly for these faults.

## Measurement

Profiling build of `rust` `abad3947e`, on the build host (i9-10900K, 20 threads, 31 GB RAM), in a window with no other heavy jobs (CPU use of every process logged with `pidstat`). Medians of 3 runs, 20 threads, `perf record` at 1234 Hz:

| Run | Wall | User | System | Peak RSS | Page faults | Kernel share of CPU | `clear_page_erms` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `ancestral --representation=dense` `rsv/a/2000` | 3.39 s | 34.4 s | 9.3 s | 8.02 GB | 2.00 million | 21.4% | 8.5% |
| `optimize --representation=dense` `rsv/a/500` | 3.80 s | 49.8 s | 5.4 s | 1.87 GB | 1.20 million | 10.5% | 3.8% |

- In the `rsv/a/2000` run, page faults (`asm_exc_page_fault`, inclusive) take 21.3 percent of CPU. 8.65 percentage points of that are under `dgemm_beta_HASWELL`: OpenBLAS sets the output array of each `cblas_dgemm` call to zero before it accumulates, and that `memset` is the first write to the new array
- 2.00 million faults of 4 KiB are 8.2 GB, the same as the peak RSS. Each page is touched for the first time once, so these faults follow from the size of the data, not from memory that is freed and touched again (for that case, see [N-alloc-jemalloc-keeps-freed-pages-per-arena.md](N-alloc-jemalloc-keeps-freed-pages-per-arena.md))
- In an earlier run on the same host while other jobs used most of its memory (5.9 GB free, 22 GB page cache), the kernel took 32.6 percent of CPU: page allocation fell into direct reclaim (`__alloc_pages_slowpath` -> `try_to_free_pages`), and the zone lock added 9.4 percent of spin-lock contention (`native_queued_spin_lock_slowpath`). In the quiet window this contention was 0.5 percent

### Transparent huge pages make the runs slower

jemalloc can mark its memory for transparent huge pages (`thp:always`; the host runs THP in `madvise` mode). Same window, 3 runs each:

| Run | jemalloc | Wall | User | System | Page faults |
| --- | --- | ---: | ---: | ---: | ---: |
| `ancestral` `rsv/a/2000` | default | 3.29-3.43 s | 34.3 s | 8.9-9.4 s | 2.00 million |
| `ancestral` `rsv/a/2000` | `thp:always` | 4.63-6.42 s | 35.6-57.2 s | 12.4-49.5 s | 32-176 thousand |
| `optimize` `rsv/a/500` | default | 3.71-3.75 s | 49.2-49.8 s | 5.3-5.5 s | 1.19-1.21 million |
| `optimize` `rsv/a/500` | `thp:always` | 4.33-4.34 s | 60.4-60.7 s | 6.5-6.7 s | 220-226 thousand |

Huge pages remove up to 98 percent of the faults, but the runs take 16 to 90 percent longer wall time and peak RSS grows by 9 to 15 percent. The kernel compaction daemon (`kcompactd0`) ran at 30 to 75 percent of one core during these runs.

## Location

- `struct DenseNodeState` [packages/treetime/src/partition/storage/dense.rs#L11-L14](../../packages/treetime/src/partition/storage/dense.rs#L11-L14) holds the node profile, `struct DenseEdgeBackward` [packages/treetime/src/partition/storage/dense.rs#L68-L71](../../packages/treetime/src/partition/storage/dense.rs#L68-L71) holds `msg_to_parent` and `msg_from_child`, and `struct DenseEdgeForward` [packages/treetime/src/partition/storage/dense.rs#L74-L76](../../packages/treetime/src/partition/storage/dense.rs#L74-L76) holds `msg_to_child`. Each is an L-by-s `Array2<f64>` in a `DenseSeqDistribution`
- For `rsv/a/2000` (3996 nodes, 3995 edges, L = 15225, s = 4), the 15981 arrays hold 7.79 GB, which accounts for almost all of the 8.02 GB peak RSS
- The products that create these arrays return a new `Array2` from `.dot()`: `fn GTR.evolve()` [packages/treetime/src/gtr/gtr.rs#L154](../../packages/treetime/src/gtr/gtr.rs#L154), `fn GTR.propagate_profile()` [packages/treetime/src/gtr/gtr.rs#L174](../../packages/treetime/src/gtr/gtr.rs#L174), and `fn in_eigenbasis()` [packages/treetime/src/gtr/gtr.rs#L259-L273](../../packages/treetime/src/gtr/gtr.rs#L259-L273)

## Memory of the large example datasets

The same count, one L-by-s array of 8-byte values per node and three per edge with s = 4, gives the dense profile memory of the large datasets:

| Dataset | Nodes | L | Profile arrays |
| --- | ---: | ---: | ---: |
| `dengue/2000` | 3906 | 10723 | 5.4 GB |
| `rsv/a/2000` | 3996 | 15225 | 7.8 GB |
| `mpox/clade-ii/100` | 198 | 197209 | 5.0 GB |
| `mpox/clade-ii/500` | 924 | 197209 | 23.3 GB |
| `mpox/clade-ii/1000` | 1916 | 197209 | 48.3 GB |
| `mpox/clade-ii/2000` | 3976 | 197209 | 100.3 GB |
| `sc2/2844` | 4455 | 29903 | 17.0 GB |
| `sc2/5100` | 10199 | 29903 | 39.0 GB |

On the build host, dense `ancestral` and `optimize` stop with an allocation failure on `mpox/clade-ii/{500,1000,2000}` and `sc2/{2844,5100}`, at 14 to 16 GB of resident memory, when the smoke test limits each process to 18.5 GiB of data segment. `dev/smoke.toml` skips these datasets in its dense rows, because whether they fit depends on the memory of the host.

## Fix direction

Reduce the number of L-by-s arrays that a run holds at the same time, for example by dropping a message after the last pass that reads it, or by computing a node profile only when an output needs it. Fewer arrays mean fewer first-touch faults and a lower peak RSS. Column compression ([N-optimize-dense-iteration-slow.md](N-optimize-dense-iteration-slow.md)) reduces L and therefore every array.

> [!IMPORTANT]
> **Investigation required.** Find which of the four arrays per node the backward pass, the forward pass, the branch-length optimization, and the outputs read, and when each is read for the last time. Then measure peak RSS, page faults, and wall time of the change on `ancestral --representation=dense rsv/a/2000` and `optimize --representation=dense rsv/a/500` at 1 and 20 threads.

## Related

- [N-io-large-dataset-memory-constraint.md](N-io-large-dataset-memory-constraint.md): the input sequences are a separate part of the peak memory
- [N-gtr-small-profile-products-pack-through-openblas.md](N-gtr-small-profile-products-pack-through-openblas.md): the other half of the time in `cblas_dgemm`
