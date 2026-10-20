# Parallel sparse leaf setup validation covers only a narrow success path

> [!IMPORTANT]
> **Investigation required.** Current evidence supports deterministic success on the exercised inputs, but correctness on multiple partitions, larger trees, and real datasets remains unverified.

## Current evidence

`fn test_marginal_sparse_parallel_pipeline_is_thread_count_deterministic()` [packages/treetime/src/ancestral/__tests__/test_marginal_sparse.rs#L499-L506](../../packages/treetime/src/ancestral/__tests__/test_marginal_sparse.rs#L499-L506) compares one-thread and four-thread results for a four-leaf tree with one nucleotide partition. It compares likelihood state and reconstructed node and edge data, so it exercises the parallel setup and indexed passes together.

A failed leaf setup leaves no partial state: `fn attach_seqs_to_graph()` takes the node inputs by value, reads the graph without mutating it, and extends the partition only after every leaf succeeds [packages/treetime/src/partition/fitch/passes.rs#L52-L87](../../packages/treetime/src/partition/fitch/passes.rs#L52-L87). In production it runs only inside `fn create_fitch_partition()`, which builds a new partition and discards it on error.

## Missing coverage

- Multiple partitions and larger leaf sets are not exercised directly.
- No automated check compares real-dataset outputs across thread counts.
- Existing dense-sparse and sparse root-invariance defects limit the strength of broader representation properties.

## Research required

Define the intended validation contract before adding tests. The contract should distinguish deterministic parallel execution from existing sparse-model divergences, then select unit, property, and golden-master coverage for each observable signal.

## Related issues

- [M-ancestral-dense-sparse-divergence.md](M-ancestral-dense-sparse-divergence.md)
- [M-ancestral-sparse-root-invariance.md](M-ancestral-sparse-root-invariance.md)
- [N-ancestral-fitch-site-classification-parallel-scaling-unverified.md](N-ancestral-fitch-site-classification-parallel-scaling-unverified.md)
