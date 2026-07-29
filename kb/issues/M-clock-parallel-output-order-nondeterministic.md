# Parallel clock CSV row order is nondeterministic

Repeated runs of the same binary record differences in `rerooted.clock.csv` and the derived SVG. `fn gather_clock_regression_results()` pushes rows into `ArrayQueue` during parallel traversal and collects queue order directly; CSV serialization and SVG rendering preserve that scheduling-dependent order [packages/treetime/src/clock/rtt.rs#L28-L80](../../packages/treetime/src/clock/rtt.rs#L28-L80).

The one- versus sixteen-thread comparisons in [kb/reports/2026-07-13_perf-report-parallel-scaling/mpox-2000-parallel-sparse-leaf-setup.md](../reports/2026-07-13_perf-report-parallel-scaling/mpox-2000-parallel-sparse-leaf-setup.md#output-equivalence) and [kb/reports/2026-07-13_perf-report-parallel-scaling/mpox-2000-zstd.md](../reports/2026-07-13_perf-report-parallel-scaling/mpox-2000-zstd.md#output-equivalence) reproduce the row-order defect. A later comparison sorted the rows and found identical fitted values across worker counts, tracing both output differences to ordering [kb/reports/2026-07-29_perf-report-work-first-dag-traversal/report.md](../reports/2026-07-29_perf-report-work-first-dag-traversal/report.md).

## Potential solutions

- O1. Produce keyed rows in parallel and sort by the approved topology/reference order before serialization.
- O2. Preallocate indexed result slots and fill each node’s deterministic slot. This avoids a sort when the ordering index already exists.

## Recommendation

Use the existing topology/reference-order index to construct deterministic slots. The row order must match the command's selected topology order and remain stable across worker counts.

## Related issues

- [N-reroot-missing-min-dev-end-to-end-oracle.md](N-reroot-missing-min-dev-end-to-end-oracle.md)
