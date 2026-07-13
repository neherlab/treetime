# iter_children_arc scans full node list per visited node

`packages/treetime-graph/src/graph_traverse.rs` `iter_children_arc` iterates the entire `self.nodes` vector for each visited node to find children, making `iter_breadth_first_backward` O(V²). Pre-existing code, not introduced by the reroot branch. Fix: traverse via outbound edge keys from each node (O(V) total). Affects `compute_div_stats` and any other caller of the serial backward traversal on large trees.
