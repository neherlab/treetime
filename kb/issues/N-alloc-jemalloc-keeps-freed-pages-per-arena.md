# jemalloc keeps freed pages per arena, so peak RSS overstates live memory

The Linux builds of the CLI and the API server allocate through jemalloc with its default configuration. jemalloc keeps freed pages resident for its decay time, separately in each arena. When the work of a run spreads over all worker threads, the freed pages of every thread's arena stay resident at the same time, and peak resident memory (peak RSS) is much higher than the memory that the run actually holds. A memory limit (cgroup, container, or batch-job limit) counts these pages like live data.

## Allocator setup

- `static GLOBAL` in [packages/app-cli/src/bin/treetime.rs#L1-L8](../../packages/app-cli/src/bin/treetime.rs#L1-L8) and in [packages/app-server/src/bin/treetime-server.rs#L1-L8](../../packages/app-server/src/bin/treetime-server.rs#L1-L8): `tikv_jemallocator::Jemalloc` as the global allocator on x86_64 and aarch64 Linux (gnu and musl). Other targets and the N-API module (`packages/app-napi`) use the system allocator
- No `malloc_conf` is set: the binaries export no `_rjem_malloc_conf` symbol, and no build sets `JEMALLOC_SYS_WITH_MALLOC_CONF`. Only `_RJEM_MALLOC_CONF` at run time changes the configuration (`docs/dev/developer_guide.md`, section "Heap profiles")
- jemalloc 5.3 defaults (`tikv-jemalloc-sys` 0.6.1): `dirty_decay_ms` 10000 (freed pages stay resident for about 10 s), `muzzy_decay_ms` 0, and 4 automatic arenas per CPU (`narenas_ratio` 4). jemalloc assigns each thread to one arena, and an arena reuses only its own freed pages

## Measurement

`mugration` on a caterpillar tree with 100,000 leaves (each internal node has one leaf child and one internal child) and one attribute with 3 states, at 20 worker threads, peak RSS in MB (default configuration: median and range of 6 runs; other rows: 2 runs):

| `_RJEM_MALLOC_CONF` | `602a051bb` (graph pass) | Tree scans |
| --- | --- | --- |
| default | 839 [827-897] | 1221 [1129-1301] |
| `dirty_decay_ms:0,muzzy_decay_ms:0` | 734 | 532 |
| `narenas:1` | 894 | 660 |

- With freed pages returned at once, or with one arena, the tree scans need 26 to 28 percent less memory than the graph pass. Only the default configuration reverses the order
- Over the GTR refinement loop, resident memory stays between 550 and 780 MB with the graph pass and grows from about 600 to 1150-1180 MB with the tree scans (20 ms samples)
- Page-fault totals differ by 2 percent (1.31 and 1.34 million), which fits retained pages and not more data
- Inferred cause: the leaf tasks of `fn GraphSchedule.scan_up()` [packages/treetime-graph/src/schedule_scan.rs#L42-L68](../../packages/treetime-graph/src/schedule_scan.rs#L42-L68) run on every worker thread, so allocations and frees spread over many arenas, and each arena keeps its freed pages for the decay time. The tree scans are also 32 percent faster on this case (wall time 10.95 to 7.47 s)

### Timetree validation runs

On the 36 single-threaded timetree runs of the simulated validation subset, release builds of the code with the forward-posterior trim against `rust` `d86d3eb65` (with the tree scans): 12 runs have more than 1.1 times the peak RSS of `rust`, at most 1.193 (42 to 50 MB), all small coalescent runs of 42 to 60 MB that grow by 4 to 10 MB, while the live heap at the peak of the worst run is 38.1 MB against 35.0 MB. The largest run of the subset grows from 123 to 133 MB (1.084). The maintainers accepted this residual for the trimming change ([kb/decisions/timetree-forward-posterior-trimming.md](../decisions/timetree-forward-posterior-trimming.md#results)); the allocator setting below remains open here.

- `_RJEM_MALLOC_CONF=lg_extent_max_active_fit:8`: all 36 runs at most 1.098 times the peak RSS of `rust`, at unchanged CPU time
- `_RJEM_MALLOC_CONF=dirty_decay_ms:200`: helps short runs only

## Fix direction

Set the decay or arena configuration of jemalloc in the binaries, for example a short `dirty_decay_ms` or fewer arenas (`narenas`). Embed the setting at build time with the `_rjem_malloc_conf` symbol or `JEMALLOC_SYS_WITH_MALLOC_CONF`, so every run uses it without an environment variable.

> [!IMPORTANT]
> **Investigation required.** The time cost of these settings is not measured. Returning freed pages at once adds `madvise` calls and page faults when the memory is used again, and fewer arenas add lock contention between worker threads. Measure wall time, user time, and peak RSS for each candidate setting on the caterpillar `mugration` case and on large standard datasets (for example `timetree dengue/2000`, `ancestral mpox/clade-ii/2000`, `ancestral --dense=true rsv/a/2000`) at 1 and 20 threads, and keep a setting only when it does not slow these runs down.
