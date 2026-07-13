# Sparse reconstruction still uses remove/insert pattern

Sparse sequence reconstruction in `marginal_sparse.rs` uses `partition.nodes.remove()` + `partition.nodes.insert()`. If an operation between removal and reinsertion panics, the entry is permanently lost from the partition map.

The marginal backward and forward pass instances of this pattern were resolved by topology-indexed pass storage. Their node and edge values move into stable slots before traversal and return to the maps after traversal completes.

## Locations

- `packages/treetime/src/partition/marginal_sparse.rs` - `reconstruct_node_sequence`: node remove/insert. A simple `get_mut` does not fit because the function reads the parent entry while writing the child entry in the same map.

## Suggested fix

Compute the reconstructed child data from immutable parent and child borrows, then acquire the mutable child entry only for the final assignment.
