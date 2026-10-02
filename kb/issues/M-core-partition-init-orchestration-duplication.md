# Partition creation is hardcoded instead of configured as complete partition sets

`build_marginal_partition()` centralizes one dense/sparse/model branch, but callers still provide construction mechanics one partition at a time and independently decide attachment, index assignment, trait-object ownership, and output handling.

## Current state

- `fn build_marginal_partition()` accepts a representation, model, graph, partition index, alphabet, and node inputs, and returns one seeded `MarginalReconstruction` [`packages/treetime/src/partition/create.rs#L34`](../../packages/treetime/src/partition/create.rs#L34). Prune builds its sparse partition with `fn build_sparse_reconstruction()` [`packages/treetime/src/partition/create.rs#L68`](../../packages/treetime/src/partition/create.rs#L68).
- Ancestral, optimize, prune, and timetree pipelines each build exactly one sequence partition and hardcode its index `0` [`packages/treetime/src/ancestral/pipeline.rs#L39`](../../packages/treetime/src/ancestral/pipeline.rs#L39) [`packages/treetime/src/optimize/pipeline.rs#L61`](../../packages/treetime/src/optimize/pipeline.rs#L61) [`packages/treetime/src/prune/pipeline.rs#L47`](../../packages/treetime/src/prune/pipeline.rs#L47) [`packages/treetime/src/timetree/pipeline.rs#L304`](../../packages/treetime/src/timetree/pipeline.rs#L304).
- Multi-partition infrastructure exists, while command configuration still creates one sequence partition.

## Impact

Adding multi-segment genomes, codon partitions, or configuration-file partitions requires coordinated changes in every pipeline.

## Required boundary

A partition configuration layer must consume complete partition descriptions and own index assignment, dense/sparse selection, GTR resolution, alignment attachment, and conversion into the capabilities requested by a pipeline. Configuration sources produce descriptions; they do not duplicate construction policy. Scientific defaults and automatic partitioning remain separate decisions.

## Related issues

- [H-core-command-module-shared-ops-entanglement.md](H-core-command-module-shared-ops-entanglement.md)
- [N-optimize-multi-alignment-input.md](N-optimize-multi-alignment-input.md)
- [N-io-multi-segment-genome-input.md](N-io-multi-segment-genome-input.md)
- [N-representation-infer-dense-stub.md](N-representation-infer-dense-stub.md)
