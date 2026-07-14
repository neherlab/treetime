# Sparse reconstruction can lose node state and suppress invariant failures

Sparse sequence reconstruction removes the child entry before resolving its parent and edge. Missing topology or partition data returns `None` while the child is absent, so the entry is permanently lost. The traversal caller then converts `None` to an empty sequence and continues successfully, hiding the invariant failure from the command.

The backward and forward passes use topology-indexed storage and do not share this failure mode.

## Evidence

- `fn PartitionMarginalSparse::reconstruct_node_sequence()` [packages/treetime/src/partition/marginal_sparse.rs#L516](../../packages/treetime/src/partition/marginal_sparse.rs#L516) removes `self.nodes[node.key]` before the fallible parent, edge, and parent-node lookups and returns `Option<Seq>`.
- `fn reconstruct_ancestral_sequences()` [packages/treetime/src/ancestral/marginal.rs#L41](../../packages/treetime/src/ancestral/marginal.rs#L41) replaces a missing reconstruction with `seq![]` and only logs a warning.
- `trait PartitionMarginalOps` [packages/treetime/src/partition/traits.rs#L163](../../packages/treetime/src/partition/traits.rs#L163) cannot distinguish intentional leaf omission from a missing node or broken topology because both are `None`.

## Required contract

- Return `Result<Option<Seq>, Report>` from the partition reconstruction API. `Ok(None)` means only that leaf output was intentionally excluded; missing nodes, parents, or edges are errors.
- Resolve every fallible lookup and compute the reconstructed sequence from immutable data before mutating the child entry.
- Commit the sequence and composition together through one final mutable child lookup. An error leaves the complete partition unchanged.
- Propagate reconstruction errors through the graph traversal. Never substitute an empty sequence for a failed reconstruction.

## Recommendation

Implement the required contract directly. It makes invalid state observable, preserves node-map completeness on every error path, and keeps the leaf-omission case explicit.
