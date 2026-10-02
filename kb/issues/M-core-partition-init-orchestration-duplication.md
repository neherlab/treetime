# Partition creation is hardcoded instead of configured as complete partition sets

`build_marginal_partition()` centralizes one dense/sparse/model branch, but callers still provide construction mechanics one partition at a time and independently decide attachment, index assignment, trait-object ownership, and output handling.

## Current state

- `fn build_marginal_partition()` accepts a representation, model, graph, partition index, alphabet, and node inputs [`packages/treetime/src/partition/create.rs#L19`](../../packages/treetime/src/partition/create.rs#L19).
- Ancestral, optimize, prune, and timetree pipelines convert the result into their own partition collections and hardcode sequence partition index `0` [`packages/treetime/src/ancestral/pipeline.rs#L132`](../../packages/treetime/src/ancestral/pipeline.rs#L132) [`packages/treetime/src/optimize/pipeline.rs#L71`](../../packages/treetime/src/optimize/pipeline.rs#L71) [`packages/treetime/src/prune/pipeline.rs#L38`](../../packages/treetime/src/prune/pipeline.rs#L38) [`packages/treetime/src/timetree/pipeline.rs#L458`](../../packages/treetime/src/timetree/pipeline.rs#L458).
- `fn initialize_partitions()` retains a second public three-way constructor in the timetree command module, but has no production caller [`packages/treetime/src/commands/timetree/initialization.rs#L99`](../../packages/treetime/src/commands/timetree/initialization.rs#L99).
- Multi-partition infrastructure exists, while command configuration still creates one sequence partition.

## Impact

Adding multi-segment genomes, codon partitions, or configuration-file partitions requires coordinated changes in every pipeline. The stale initializer can also drift from the active creator without compiler or test pressure.

## Required boundary

A partition configuration layer must consume complete partition descriptions and own index assignment, dense/sparse selection, GTR resolution, alignment attachment, and conversion into the capabilities requested by a pipeline. Configuration sources produce descriptions; they do not duplicate construction policy. Scientific defaults and automatic partitioning remain separate decisions.

Delete the unused `initialize_partitions()` path as part of consolidation.

## Related issues

- [H-core-command-module-shared-ops-entanglement.md](H-core-command-module-shared-ops-entanglement.md)
- [N-optimize-multi-alignment-input.md](N-optimize-multi-alignment-input.md)
- [N-io-multi-segment-genome-input.md](N-io-multi-segment-genome-input.md)
- [N-representation-infer-dense-stub.md](N-representation-infer-dense-stub.md)
