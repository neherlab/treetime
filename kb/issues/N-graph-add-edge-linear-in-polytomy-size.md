# `Graph.add_edge()` is linear in the number of children of the source

`fn Graph.add_edge()` checks every outbound edge of the source node for an existing edge to the same target before it adds the edge.

## Location

[`packages/treetime-graph/src/graph_ops.rs#L42-L60`](../../packages/treetime-graph/src/graph_ops.rs#L42-L60)

## Impact

Adding $k$ children to one node costs $O(k^2)$. Trees built from low-divergence data have polytomies with hundreds to thousands of children, and the tree builder attaches samples and moves subtrees there repeatedly. `fn Graph.remove_edge()` and `fn Graph.remove_node()` were made proportional to the degree for the same reason on branch `tree-builder`.

## Fix

Check the inbound edges of the target instead: in a tree the target has at most one parent, so the check becomes constant-time. Alternatively, move the check to debug builds.
