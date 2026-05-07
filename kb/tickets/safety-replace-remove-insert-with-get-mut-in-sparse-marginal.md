# Replace remove/insert with get_mut in sparse marginal passes

The sparse marginal backward and forward passes in `marginal_passes.rs` and `marginal_sparse.rs` use `partition.edges.remove()` + `partition.edges.insert()` (and one `partition.nodes.remove()` + `insert()`). If any operation between removal and reinsertion panics, the entry is permanently lost from the partition map.

The dense forward pass already uses shared borrows plus `get_mut` for the equivalent update, exploiting disjoint field borrows of `partition.nodes` and `partition.edges`.

## Locations

- `packages/treetime/src/representation/partition/marginal_passes.rs:107:` - `process_node_backward`: edge remove/insert. Fixable with `get_mut` (computation uses only `partition.gtr` and locals, no borrow conflict with `partition.edges`).
- `packages/treetime/src/representation/partition/marginal_passes.rs:136:` - `process_node_forward` parent loop: edges removed, cloned into backup Vec, re-inserted at lines 170-172. More complex because the code reads from multiple edge entries while mutating the collection. May require collecting read data before the mutable pass.
- `packages/treetime/src/representation/partition/marginal_passes.rs:189:` - `process_node_forward` child loop: edge remove/insert. Structurally identical to the dense path that already uses shared borrows. Disjoint field borrows (`&partition.nodes` + `&mut partition.edges`) apply directly.
- `packages/treetime/src/representation/partition/marginal_sparse.rs:382:` - `reconstruct_node_sequence`: node remove/insert. Cannot use simple `get_mut` because the code reads `self.nodes[parent_key]` (line 388) while needing to write `self.nodes[node.key]` (line 423). Two entries in the same map require either temporary extraction or computing the result into locals first.

## Suggested fix

Apply the same disjoint-field-borrow + `get_mut` pattern for locations 1 and 3. Locations 2 and 4 need restructuring to separate read and write phases.

## Related issues

- Source: [N-ancestral-sparse-remove-insert-pattern.md](../issues/N-ancestral-sparse-remove-insert-pattern.md) -- delete after full resolution
