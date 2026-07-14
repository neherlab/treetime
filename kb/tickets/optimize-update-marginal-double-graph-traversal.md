# update_marginal traverses the graph twice for mixed partitions

`update_marginal()` is generic over `P: PartitionMarginalOps<N, E>`. Dense and sparse partitions are different concrete types, so each call requires a separate slice. Every site that needs both runs two sequential graph traversals:

- `packages/treetime/src/commands/optimize/run.rs:173+180:` initial marginal
- `packages/treetime/src/commands/optimize/run.rs:221+222:` pre-annotation
- `packages/treetime/src/commands/optimize/run.rs:314+315:` numerical failure revert
- `packages/treetime/src/commands/optimize/run.rs:352+353:` worsened revert
- `packages/treetime/src/commands/optimize/run.rs:395+396:` per-iteration likelihood

Each traversal does a full backward + forward pass over the tree. The per-node partition work dominates, but the graph traversal overhead and synchronization barriers double.

A single-traversal version could accept both partition types via a trait object (`dyn PartitionMarginalOps`) or an enum wrapper, processing all partitions at each node in one pass.

## Related issues

- Source: [kb/issues/N-optimize-double-graph-traversal-update-marginal.md](../issues/N-optimize-double-graph-traversal-update-marginal.md) -- delete after full resolution
