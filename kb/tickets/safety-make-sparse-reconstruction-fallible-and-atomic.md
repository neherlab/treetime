# Make sparse reconstruction fallible and atomic

Make ancestral reconstruction report invariant failures and preserve partition state when reconstruction cannot complete.

## Required changes

- Change `PartitionMarginalOps::reconstruct_node_sequence()` to return `Result<Option<Seq>, Report>` in [packages/treetime/src/partition/traits.rs#L163](../../packages/treetime/src/partition/traits.rs#L163). Reserve `Ok(None)` for a leaf intentionally excluded by `include_leaves`.
- Update the dense implementation to report a missing node as an error instead of returning `None`.
- In `PartitionMarginalSparse::reconstruct_node_sequence()` [packages/treetime/src/partition/marginal_sparse.rs#L516](../../packages/treetime/src/partition/marginal_sparse.rs#L516), resolve the child, parent key, edge key, parent data, and edge data before any mutation.
- Compute the reconstructed sequence and composition in local values, then update the child through one final `get_mut()` after every fallible operation has succeeded.
- In `reconstruct_ancestral_sequences()` [packages/treetime/src/ancestral/marginal.rs#L41](../../packages/treetime/src/ancestral/marginal.rs#L41), propagate the error with node context. Remove the empty-sequence fallback and warning.

## Validation

- Missing child, parent, edge, and parent-partition data each return an error containing the node key and leave `nodes` byte-for-byte unchanged.
- Intentional leaf exclusion returns `Ok(None)` without mutation.
- Successful root, internal-node, and included-leaf reconstruction returns `Ok(Some(sequence))` and updates sequence plus composition together.
- The traversal returns the reconstruction error and does not call the visitor with an empty sequence.
- Run the full lint and test suite.

## Related issues

- Source: [kb/issues/N-ancestral-sparse-remove-insert-pattern.md](../issues/N-ancestral-sparse-remove-insert-pattern.md)
