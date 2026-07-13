# Profiling: Why Sparse-Mode Passes Scale Poorly

**Commit:** `3f1820a5ee5f443691624e00ab4aad4e4012a049`
**Method:** `perf record --call-graph dwarf` on the `profiling` build (release + debug symbols), run on host. Dataset `data/mpox/clade-ii/2000` (~2000 tips, 197 kb alignment).

## Summary

The modest parallel speedups are not caused by the parallel passes being slow. They are caused by **how little of the total work is actually inside those passes**, and by **fork-join barrier overhead** from parallelising one tree depth level at a time.

Two distinct signatures:

- **Ancestral** (1.26x at 16 threads): dominated by serial work outside the parallel passes. Only ~26% of runtime is parallel compute.
- **Mugration** (2x at 16 threads, but 6+ cores busy): real per-node compute, but ~21% of thread time is spent spinning in the work-stealing scheduler between depth levels.

## Root cause: per-level fork-join

Both sparse ancestral and mugration traverse the tree one BFS depth level ("frontier") at a time. The pass driver (`packages/treetime/src/partition/indexed_pass.rs:184`) is:

```rust
for range in self.frontiers.clone() {   // serial loop over depth levels
    let frontier = &mut slots[range];
    visit(..., frontier)?;               // calls frontier.par_iter_mut() inside
}
```

Every depth level is a separate `par_iter_mut()` with an implicit rayon barrier at its end. The mpox-2000 tree has **121 depth levels** (mpox-500: 38), average frontier width 33 nodes (mpox-500: 24). So per pass, rayon forks and joins ~121 times, each time over a few dozen nodes.

Item count (24-33) is enough to hand every thread something. The problem is the **ratio of useful work per level to fork-join cost per level**: with narrow levels and, in sparse mode, small per-node work, threads finish their slots almost immediately and then spin or sleep at the barrier waiting for the level to close.

## Ancestral: serial work dominates

Single-threaded self-time (`tmp/profiling/anc-j1.data`, 7557 samples), categorised:

| Category                        | Share | Notable symbols                                                                   |
| ------------------------------- | ----: | --------------------------------------------------------------------------------- |
| Parallel compute (fitch passes) | 25.6% | `resolve_fixed_positions_backward` (22%), `run_fitch_forward_indexed`             |
| Serial input I/O                | 18.8% | `FastaReader::read`, `lzma_crc64`, `lzma_decode`                                  |
| Serial setup                    | 17.2% | `SparseNodePartition::new`                                                        |
| Serial output                   |  7.5% | `serde_json format_escaped_str`, `write_augur_node_data_json`, `BufWriter::flush` |
| Kernel page mgmt + memcg + misc | 30.6% | page faults, `clear_page`, `memmove`, orchestration                               |

The serial setup is `attach_seqs_to_graph` (`packages/treetime/src/ancestral/fitch.rs:55`), a plain serial loop over leaves that builds each leaf's sparse partition (`SparseNodePartition::new`). It is 17% of runtime and is embarrassingly parallel (one independent sparse representation per leaf), but runs on one thread.

Even with a perfectly scaling parallel section, Amdahl's law with a serial fraction of ~0.5 caps speedup at ~2x. Measured 1.26x, because the 26% parallel part also pays the per-level barrier cost above.

The j=8 profile (`anc-j8.data`) is nearly identical in relative shares: the serial functions do not shrink, and rayon overhead stays small (`bridge_producer_consumer` 0.19%) because there is so little parallel work to fight over.

## Allocator: jemalloc is in use and is not a bottleneck

The binary links jemalloc as the global allocator (`packages/app-cli/src/bin/treetime.rs:7`, gated to Linux gnu/musl on x86_64/aarch64, active in the profile: `_rjem_malloc`, `do_rallocx`). Breaking down the "kernel page mgmt + memcg + misc" bucket from the ancestral j=1 profile:

| Sub-category                             | Share | Symbols                                                                           |
| ---------------------------------------- | ----: | --------------------------------------------------------------------------------- |
| jemalloc userspace                       |  0.5% | `_rjem_malloc`, `do_rallocx`, `_rjem_je_*`                                        |
| Kernel page fault + zeroing              | ~6.0% | `clear_page_erms`, `__alloc_pages`, `do_anonymous_page`, `get_page_from_freelist` |
| cgroup memory accounting                 | ~2.1% | memcg hooks from the host's cgroup v2 memory controller (systemd)                 |
| BTreeMap ops                             |  0.4% | node inserts                                                                      |
| Other kernel + `memmove` + orchestration |  ~21% | I/O syscalls, JSON/FASTA copies, `pipeline::run`                                  |

The allocator itself is 0.5%: swapping it (mimalloc, etc.) is not worthwhile. The allocator-adjacent cost is the ~6% kernel page-fault + zeroing, driven by allocation churn: `SparseNodePartition::new` builds many small per-node `BTreeMap`s and `Vec`s and frees them, so the growing heap keeps faulting and zeroing fresh pages. Levers, by expected value:

- **L1. Cut allocation churn in sparse construction.** Reuse buffers/arenas across nodes, pre-size collections, or replace per-node `BTreeMap` with a flat sorted `Vec<(pos, val)>` where ordered-map lookup is not needed. Shrinks both the 17% construction compute and its page-fault tail at the source.
- **L2. Tune jemalloc page retention.** Raise `dirty_decay_ms`/`muzzy_decay_ms` or enable background threads via `MALLOC_CONF` so freed pages stay mapped instead of being returned to the kernel and re-faulted. Targets the ~6% page-fault cost; cheap to try.
- **L3. The ~2% memcg is host cgroup v2 overhead** (systemd puts every process in a memory cgroup; the profile ran on the host, not in a container). It tracks the same page-fault traffic as L1/L2, so cutting allocation churn reduces it too. Removing it outright needs `cgroup_disable=memory` at the host level, which is a deployment choice, not a code change.

## Mugration: barrier spin dominates the overhead

Mugration does more real per-node work: exponentiating the country-state GTR transition matrix via OpenBLAS `dgemm`. Single-threaded-per-call BLAS (confirmed `SINGLE_THREADED` / `Sequential` via `treetime debug`), so there is **no BLAS-under-rayon oversubscription**.

j=16 self-time (`mug-j16.data`, 51934 samples):

| Category                          |     Share | Symbols                                                                                                              |
| --------------------------------- | --------: | -------------------------------------------------------------------------------------------------------------------- |
| Matrix exponential (real compute) |      ~28% | `dgemm_kernel_HASWELL`, `dgemm_itcopy`, `dgemm_otcopy`                                                               |
| Work-stealing / rayon overhead    |      ~12% | `crossbeam_deque::Stealer::steal` (6.65%), `rayon bridge_producer_consumer` (3.76%), `WorkerThread::wait_until_cold` |
| Kernel scheduler churn            |     ~9.5% | `try_to_wake_up`, `native_queued_spin_lock_slowpath`, `__schedule`, `pick_next_task_fair`                            |
| Mutex                             |       ~2% | `pthread_mutex_lock`/`unlock`                                                                                        |
| ndarray + math + misc             | remainder | `Zip::map_collect_owned`, `__ieee754_exp`, `count_transitions`                                                       |

`Stealer::steal` at 6.65% is workers spinning to steal work that is not there; `try_to_wake_up` and the queued spin lock are the kernel waking and re-sleeping threads. This ~21% combined is pure parallelism machinery, paid because each of the 121 frontier levels is its own fork-join over a few dozen nodes. It explains the earlier observation of 6+ cores "busy" at only 1.3-2x speedup: much of that CPU time is spin-waiting, which `perf` counts as busy.

## Recommendations

Ordered by expected impact on sparse-mode scaling.

- **R1. Parallelise sequence-compression setup.** `attach_seqs_to_graph` (`fitch.rs:55`) is a serial 17% of ancestral. Build the per-leaf `SparseNodePartition` values in parallel (`par_iter` over leaves, collect into the map). Removes a large chunk of the serial fraction directly.
- **R2. Reduce fork-join count.** The per-level barrier is the shared root cause. Options: process independent subtrees (whose depth ranges do not interact) as one larger parallel region instead of one-level-at-a-time; or coarsen granularity so each rayon task covers a subtree, not a single node. Fewer, larger parallel regions means the barrier cost is amortised over more work.
- **R3. Overlap I/O with compute.** Ancestral spends ~19% reading and decompressing the alignment and ~7.5% writing JSON, both serial and both currently outside any parallel region. Streaming the FASTA/xz read and pipelining output serialisation would shrink the serial fraction.
- **R4. Parallelise across partitions** when multiple partitions exist (e.g. amino-acid partitions): `marginal_backward`/`marginal_forward` iterate partitions serially. No effect on the current single-nucleotide-partition runs, but relevant as soon as AA partitions land.

R1 and R3 attack ancestral's serial fraction (the binding constraint there). R2 attacks the barrier overhead that limits every frontier-parallel pass, mugration included.

## Reproduction

```bash
# build the profiling binary (release + debug symbols) in the container
./dev/docker/run ./dev/dev bp treetime

# record on host (kernel perf access required)
perf record -F 1234 --call-graph dwarf -o tmp/profiling/anc-j1.data -- \
  .build/docker/profiling/treetime ancestral --method-anc=marginal --dense=false \
  --tree=data/mpox/clade-ii/2000/tree.nwk --aln=data/mpox/clade-ii/2000/aln.fasta.xz \
  --output-all=tmp/profiling/anc-out -j 1 --silent

perf report -i tmp/profiling/anc-j1.data --stdio --no-children -g none
```
