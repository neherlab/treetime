# Parallel sparse leaf setup may regress single-thread ancestral runtime

> [!IMPORTANT]
> **More performance research is required.** The available measurement covers one dataset and does not establish whether the regression is systematic or noise.

## Observation

The mpox-2000 comparison reports ancestral wall time increasing from 5.732 seconds to 6.057 seconds at `-j 1`, a 5.7% regression. The same change improves ancestral wall time at 2-16 threads. Optimize and timetree do not show a comparable one-thread regression in that run.

The benchmark used one warmup and three measured runs per binary and configuration. That supports the recorded observation but is insufficient to choose a sequential fast path or accept a general one-thread tradeoff.

## Research required

- Repeat paired measurements across representative datasets and compression ratios.
- Separate Rayon collection overhead from changed allocation and map-construction behavior.
- Establish whether the regression exceeds normal run-to-run variation and affects end-to-end commands.
- Compare any proposed one-thread path for output and error-state equivalence before implementation.

## Evidence

- [kb/reports/2026-07-13_perf-report-parallel-scaling/mpox-2000-parallel-sparse-leaf-setup.md](../reports/2026-07-13_perf-report-parallel-scaling/mpox-2000-parallel-sparse-leaf-setup.md)

## Related issues

- [N-ancestral-parallel-sparse-leaf-validation-coverage.md](N-ancestral-parallel-sparse-leaf-validation-coverage.md)
- [N-ancestral-parallel-sparse-leaf-error-atomicity-unverified.md](N-ancestral-parallel-sparse-leaf-error-atomicity-unverified.md)
