# Parallel sparse leaf setup validation covers only a narrow success path

> [!IMPORTANT]
> **Investigation required.** Current evidence supports deterministic success on the exercised inputs, but correctness on multiple partitions, larger trees, and real datasets remains unverified.

## Current evidence

`fn test_marginal_sparse_parallel_pipeline_is_thread_count_deterministic()` [packages/treetime/src/ancestral/__tests__/test_marginal_sparse.rs#L474-L480](../../packages/treetime/src/ancestral/__tests__/test_marginal_sparse.rs#L474-L480) compares one-thread and four-thread results for a four-leaf tree with one nucleotide partition (6 sites). It compares likelihood state and reconstructed node and edge data, so it exercises the parallel setup and indexed passes together.

`fn test_fitch_pass_determinism_results_are_equal_across_thread_counts_and_runs()` [packages/treetime/src/ancestral/__tests__/test_fitch_pass_determinism.rs#L9-L27](../../packages/treetime/src/ancestral/__tests__/test_fitch_pass_determinism.rs#L9-L27) compares the Fitch compression (`fn compress_sequences()`) across one and four threads and repeated runs. Its input `fn fitch_input()` [packages/treetime/src/ancestral/__tests__/test_fitch_pass_support.rs#L29](../../packages/treetime/src/ancestral/__tests__/test_fitch_pass_support.rs#L29) is the same four-leaf tree with one 16-site nucleotide alignment.

A failed leaf setup leaves no partial state: `fn attach_seqs_to_graph()` takes the partition and the node inputs by value, reads the graph without mutating it, and extends the partition only after every leaf succeeds [packages/treetime/src/partition/fitch/passes.rs#L54-L87](../../packages/treetime/src/partition/fitch/passes.rs#L54-L87). In production it runs only inside `fn compress_sequences()`, which only `fn create_fitch_partition()` calls [packages/treetime/src/partition/fitch/passes.rs#L28-L52](../../packages/treetime/src/partition/fitch/passes.rs#L28-L52); that function builds a new partition and discards it on error.

## Missing coverage

- Multiple partitions and larger leaf sets are not exercised directly: both thread-count tests use the same four-leaf tree and one partition.
- No automated check compares real-dataset outputs across thread counts.
- Existing dense-sparse and sparse root-invariance defects limit the strength of broader representation properties.

## Research required

Define the intended validation contract before adding tests. The contract should distinguish deterministic parallel execution from existing sparse-model divergences, then select unit, property, and golden-master coverage for each observable signal.

## Related issues

- [M-ancestral-dense-sparse-divergence.md](M-ancestral-dense-sparse-divergence.md)
- [M-ancestral-sparse-root-invariance.md](M-ancestral-sparse-root-invariance.md)
- [N-ancestral-fitch-site-classification-parallel-scaling-unverified.md](N-ancestral-fitch-site-classification-parallel-scaling-unverified.md)
