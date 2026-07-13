# Parallel sparse leaf setup validation covers only a narrow success path

> [!IMPORTANT]
> **Investigation required.** Current evidence supports deterministic success on the exercised inputs, but broader correctness and failure behavior remain unverified.

## Current evidence

`fn test_marginal_sparse_parallel_pipeline_is_thread_count_deterministic()` [packages/treetime/src/ancestral/__tests__/test_marginal_sparse.rs#L544](../../packages/treetime/src/ancestral/__tests__/test_marginal_sparse.rs#L544) compares one-thread and four-thread results for a four-leaf tree with one nucleotide partition. It compares likelihood state and reconstructed node and edge data, so it exercises the parallel setup and indexed passes together.

The mpox-2000 performance report additionally records byte-identical ancestral, optimize, timetree, and mugration outputs between the baseline and changed binaries at selected thread counts. That practical check is retained as a report rather than an automated regression test.

## Missing coverage

- Multiple partitions and larger leaf sets are not exercised directly.
- Setup failures and partial mutation are not tested.
- The report's real-dataset output equivalence is not represented by a stable automated oracle.
- Existing dense-sparse and sparse root-invariance defects limit the strength of broader representation properties.

## Research required

Define the intended validation contract before adding tests. The contract should distinguish deterministic parallel execution from existing sparse-model divergences, then select unit, property, and golden-master coverage for each observable signal.

## Related issues

- [N-ancestral-parallel-sparse-leaf-error-atomicity-unverified.md](N-ancestral-parallel-sparse-leaf-error-atomicity-unverified.md)
- [M-ancestral-dense-sparse-divergence.md](M-ancestral-dense-sparse-divergence.md)
- [M-ancestral-sparse-root-invariance.md](M-ancestral-sparse-root-invariance.md)
- [N-ancestral-parallel-sparse-leaf-single-thread-regression.md](N-ancestral-parallel-sparse-leaf-single-thread-regression.md)
- [kb/reports/2026-07-13_perf-report-parallel-scaling/mpox-2000-parallel-sparse-leaf-setup.md](../reports/2026-07-13_perf-report-parallel-scaling/mpox-2000-parallel-sparse-leaf-setup.md)
