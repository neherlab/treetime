# Fitch site classification regresses parallel command performance

The indexed Fitch implementation builds `struct FitchSiteIndex` with a serial leaf-by-site scan before parallel reconstruction [[src](https://github.com/neherlab/treetime/blob/20164657c95d853f15e1a5dc7d8955c6134798ce/packages/treetime/src/ancestral/fitch.rs#L34-L73)]. Representative paired benchmarks report slower eight- and sixteen-worker ancestral, optimize, and timetree runs even though the change improves isolated single-worker work [[doc](https://github.com/neherlab/treetime/blob/20164657c95d853f15e1a5dc7d8955c6134798ce/kb/reports/2026-07-13_perf-report-parallel-scaling/mpox-2000-fitch-informative-positions.md#L9-L27)].

The issue is performance-only: tested outputs remain byte-identical and no scientific result, crash, or unavailable standard feature is implicated. It therefore uses the KB's negligible severity; the representative multi-worker impact still requires correction before accepting the indexed implementation.

## Potential solutions

- O1. Parallelize independent site classification and merge in stable site order.
- O2. Fuse classification into the existing parallel reconstruction pass. This can avoid a separate traversal but couples indexing to reconstruction state.

## Recommendation

Classify independent alignment sites in parallel, preserve deterministic site ordering, and merge per-site results once. Evaluate the complete commands across representative worker counts; a one-worker Rayon pool is not a scaling benchmark.

## Related issues

- [N-ancestral-parallel-sparse-leaf-single-thread-regression.md](N-ancestral-parallel-sparse-leaf-single-thread-regression.md)
