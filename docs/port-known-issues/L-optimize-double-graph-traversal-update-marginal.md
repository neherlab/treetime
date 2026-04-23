# update_marginal traverses the graph twice for mixed partitions

`update_marginal()` is generic over `P: PartitionMarginalOps<N, E>`. Dense and sparse partitions are different concrete types, so each call requires a separate slice. Every site that needs both runs two sequential graph traversals:

- `packages/treetime/src/commands/optimize/run.rs:174+181:` initial marginal
- `packages/treetime/src/commands/optimize/run.rs:216+217:` pre-annotation
- `packages/treetime/src/commands/optimize/run.rs:309+310:` numerical failure revert
- `packages/treetime/src/commands/optimize/run.rs:347+348:` worsened revert
- `packages/treetime/src/commands/optimize/run.rs:390+391:` per-iteration likelihood

Each traversal does a full backward + forward pass over the tree. The per-node partition work dominates, but the graph traversal overhead and synchronization barriers double.

A single-traversal version could accept both partition types via a trait object (`dyn PartitionMarginalOps`) or an enum wrapper, processing all partitions at each node in one pass.
