# optimize on graphs of segment trees keeps the topology fixed

`optimize --segment-trees --segment-alignments` fits branch lengths on a graph of segment trees, with one length per shared edge fitted to all segments that carry it ([kb/algo/optimization.md](../algo/optimization.md#branch-lengths-on-graphs-of-segment-trees)). Remaining:

- **Topology edits and rerooting** (zero-length collapse, merging of shared mutations, reversion hoisting, reroot) are off on a segment graph until the edge scopes have rules for them ([kb/proposals/edge-scoped-inheritance.md](../proposals/edge-scoped-inheritance.md#topology-operations)); a polytomy of one segment stays unresolved
- **Series chains with different scopes**: two edges in series through a node with one parent and one child keep their length ratio, as on a tree, also when they carry different segments, so the lengths of such edges are not fitted separately; the indel term of a chain uses the rate of its first edge
- **Pass-through nodes that topology does not determine** keep a segment's own edge beside a shared path, so the two lengths are fitted separately ([kb/algo/graph.md](../algo/graph.md#merging-trees-by-name))
