# Polytomy resolution duplicates topology primitives

Two code quality problems in `packages/treetime/src/timetree/optimization/polytomy.rs`:

## Duplicated find_polytomy_nodes

`find_polytomy_nodes` is implemented identically in both `polytomy.rs` and `merge_shared_mutations.rs`. Both filter graph nodes by `degree_out() > 2`. Should be a single generic function since the logic depends only on `Graph<N, E, D>` traits.

## remove_obsolete_nodes does not sum branch lengths

`remove_obsolete_nodes` collapses single-child internal nodes using `graph.collapse_edge()`, which reparents children but does not sum branch lengths into the merged edge. The correct implementation `remove_node_if_trivial` in `reroot.rs` uses `graph.remove_node()` + `graph.add_edge()` with summed branch lengths. The bug is currently masked because `prepare_tree_after_topology_change` wipes all cached distributions afterward, but branch length information is lost.
