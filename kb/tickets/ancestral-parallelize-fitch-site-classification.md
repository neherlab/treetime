# Parallelize Fitch site classification

Remove the serial `FitchSiteIndex` construction bottleneck while retaining deterministic reconstruction.

## Required changes

- Compute per-site classification independently under Rayon.
- Merge classifications in stable site order.
- Avoid shared mutable maps in workers.
- Keep binary and multifurcation scientific behavior unchanged; the separate polytomy issue governs algorithm semantics.

## Validation

- Exact one-worker/multi-worker output equivalence.
- Representative ancestral, optimize, and timetree benchmarks at one, two, eight, and sixteen workers.
- Regression threshold based on paired before/after observations grouped by revision.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/N-ancestral-fitch-site-classification-parallel-regression.md](../issues/N-ancestral-fitch-site-classification-parallel-regression.md)
