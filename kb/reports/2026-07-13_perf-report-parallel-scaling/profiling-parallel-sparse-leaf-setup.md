# Profiling Parallel Sparse Leaf Setup

**Verdict:** the intended setup work is now distributed across Rayon workers. `SparseNodePartition::new` remains a large share of aggregate CPU samples because the same construction work still exists, but it is no longer a serial wall-time component at multi-thread settings. The next ancestral compute target is the dense fixed-position backward scan.

- **Commit:** `ac1632682aa0e3c8444dc173d0112ded32491e91`
- **Dataset:** `data/mpox/clade-ii/2000`
- **Method:** `perf record -F 1234 --call-graph dwarf` on a profiling build, with the process pinned to CPUs 2–17

## Finding: setup parallelization is active

At `-j 8`, `SparseNodePartition::new` accounts for 16.97% of aggregate self samples. Its samples are balanced across all eight Rayon worker threads:

| Worker TID | Self share |
| ---------: | ---------: |
|     129254 |      2.11% |
|     129255 |      2.09% |
|     129256 |      2.05% |
|     129257 |      2.10% |
|     129258 |      2.17% |
|     129259 |      2.16% |
|     129260 |      2.17% |
|     129261 |      2.11% |

The implementation first resolves leaf records in parallel, then constructs a partition-local map with `par_iter`, and finally acquires the partition write lock once to extend the node map. See [`packages/treetime/src/ancestral/fitch.rs`](../../../packages/treetime/src/ancestral/fitch.rs#L55).

The one-thread profile still attributes 18.04% to `SparseNodePartition::new`, close to the 17.2% reported before the change in [profiling-findings.md](profiling-findings.md). This is expected: parallelization changes when the work completes, not the CPU instructions required to build each sparse node. At eight threads, the aggregate 16.97% is concurrent worker time and must not be interpreted as a 16.97% serial fraction.

The wall-time measurements corroborate the thread attribution. Ancestral improves by 20.1% at 8 threads and 21.2% at 16 threads, while CPU utilization rises from 1.38 to 1.87 at 8 threads and from 1.66 to 2.47 at 16 threads. Full measurements are in [mpox-2000-parallel-sparse-leaf-setup.md](mpox-2000-parallel-sparse-leaf-setup.md).

## Current ancestral profile

| Symbol or activity                                              | `-j 1` | `-j 8` | Interpretation                                  |
| --------------------------------------------------------------- | ------: | ------: | ----------------------------------------------- |
| `resolve_fixed_positions_backward`                              |  18.24% |  18.45% | Parallel dense alignment scan                   |
| `SparseNodePartition::new`                                      |  18.04% |  16.97% | Parallel at `-j 8`; balanced over eight workers |
| `FastaReader::read`                                             |  14.16% |  12.41% | Main-thread input                               |
| `lzma_crc64` + `lzma_decode`                                    |   5.05% |   4.60% | Main-thread xz decompression                    |
| `__memmove_avx_unaligned_erms`                                  |   7.65% |   9.69% | Cross-cutting memory movement                   |
| `write_augur_node_data_json_with_aa` + `format_escaped_str`     |   9.37% |   8.63% | Serial output formatting                        |
| `run_fitch_forward_indexed`                                     |   2.34% |   2.24% | Parallel forward pass                           |

The captures contain 7,058 samples at `-j 1` and 7,969 samples at `-j 8`, with zero lost samples. Percentages are self samples across all threads, so overlapping functions are not double-counted in the table.

## What to try next

### Ancestral: reduce the fixed-position work domain

`resolve_fixed_positions_backward` is now the largest named compute symbol at 18.45%. It scans every child sequence position for every internal node even though variable Fitch states are sparse. See [`packages/treetime/src/ancestral/fitch_sub.rs`](../../../packages/treetime/src/ancestral/fitch_sub.rs#L83).

The next ancestral experiment should compress repeated invariant alignment patterns while retaining one entry per genuinely variable position, then run substitution Fitch on that compressed domain. A lower-risk kernel comparison can first seed a parent sequence from one child and scan only the remaining children. Both target useful work rather than adding scheduling around an already parallel kernel.

### Shared graph passes: retain frontier coarsening as a separate target

The original mugration profile found substantial work-stealing and scheduler churn from one fork-join per tree frontier. That mechanism remains in [`packages/treetime/src/partition/indexed_pass.rs`](../../../packages/treetime/src/partition/indexed_pass.rs#L185). This follow-up profile covers ancestral only and therefore does not supersede the mugration evidence in [profiling-findings.md](profiling-findings.md). Dependency-aware subtree tasks or another coarsened scheduling unit remain the next shared traversal experiment.

### Serial boundary work

FASTA reading plus xz decoding accounts for about 17.0% of `-j 8` self samples on the main thread, and output formatting accounts for 8.6%. The completed zstd experiment in [mpox-2000-zstd.md](mpox-2000-zstd.md) found no measurable end-to-end one-thread gain from codec replacement alone, so the remaining I/O hypothesis is overlap or reduced data movement rather than another compression-format substitution.

## Reproduction

```bash
./dev/docker/run taskset -c 2-17 ./dev/dev bp treetime

taskset -c 2-17 perf record -F 1234 --call-graph dwarf \
  -o tmp/parallel-sparse-leaf-setup/profiling/ancestral-j8.data -- \
  .build/docker/profiling/treetime ancestral \
  --method-anc=marginal --dense=false \
  --tree=data/mpox/clade-ii/2000/tree.nwk \
  --aln=data/mpox/clade-ii/2000/aln.fasta.xz \
  --output-all=tmp/parallel-sparse-leaf-setup/profiling/output-j8 \
  -j 8 --silent

perf report \
  -i tmp/parallel-sparse-leaf-setup/profiling/ancestral-j8.data \
  --stdio --no-children -g none
```
