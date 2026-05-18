# Extract prune domain logic from commands

## Description

`commands/prune/run.rs` contains tree-pruning algorithms that belong in a domain module. No `src/prune/` module exists.

Domain functions in `commands/prune/run.rs`:

- `prune_nodes` (pub(super), 12 lines): orchestrates internal-node pruning then leaf pruning
- `prune_internal_nodes` (37 lines): identifies edges to collapse by short branch length, zero mutations, or name match
- `prune_leaves` (25 lines): identifies leaf edges to collapse by name match
- `collapse_sparse_edges_from_leaf_recursive` (pub(super), 23 lines): collapses leaf edge and propagates upward through single-child ancestors
- `get_edge_num_muts` (pub(super), 16 lines): counts mutations across partitions for one edge
- `should_collapse_parent` (3 lines): checks single-child non-root condition

These are tree surgery operations. They take graph + partitions, not CLI args.

Create `src/prune/` module with these functions. `commands/prune/run.rs` becomes a thin orchestrator: parse args, load tree, call `prune::prune_nodes()`, call `merge_shared_mutation_branches()`, write output.

## Validation

- `src/prune/` compiles independently with no `commands/` imports
- All existing tests pass (including `commands/prune/__tests__/`)
- `commands/prune/run.rs` imports from `crate::prune` instead of defining functions locally

## Related issues

Source: [H-core-multi-client-architecture-library-purity](../issues/H-core-multi-client-architecture-library-purity.md)
