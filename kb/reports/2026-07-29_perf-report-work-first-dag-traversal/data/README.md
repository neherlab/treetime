# Raw data

- `binary-hashes.txt`: SHA-256 of the four binaries (before/after x release/profiling).
- `hyperfine/*.json`: wall-clock measurements, one file per (command, dataset). Each result carries its parameters (`jobs`, and `rev` or the full binary path under `bin`) and the per-run `times` array, so revision/dataset/jobs/replicate are all recoverable.
- `perfstat/*.csv`: `perf stat -x,` performance counters at thread endpoints {1, 16} for timetree, sparse marginal, and parsimony on dengue/2000, both revisions.
- `timev/*.txt`: `/usr/bin/time -v` resource output for the same points.
- `profile/*.top-symbols.txt`: `perf report` flat and call-graph top symbols for the four profiled points (timetree j8 before/after, sparse marginal j8 after, parsimony j16 after). The raw `perf.data` files (hundreds of MB) are not committed; regenerate with `../scripts/perf_capture.sh` while the system is idle.

## Scripts

- `run_one.sh` / `run_bench.sh`: command dispatch (command + dataset + jobs -> treetime invocation; revision label -> binary).
- `bench_hyperfine.sh`, `bench_hyperfine_dense.sh`: the sparse and dense wall-clock matrices.
- `perf_capture.sh`: `perf stat` + `perf record` capture set (run while the system is idle).
- `analyze.py`: aggregates hyperfine JSON (and perf stat CSV) into runtime, scaling, and performance-counter tables with bootstrap CIs.
- `plot_report.py`: renders the three report charts (SVG + PNG) into `../assets/` from the hyperfine JSON.
