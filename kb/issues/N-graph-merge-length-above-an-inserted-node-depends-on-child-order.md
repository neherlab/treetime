# The merge takes the length above an inserted node from one child of a polytomy

When several children of a polytomy `X` of segment B move below the node `W` that another segment puts above them, the first child that the merge visits creates the chain `X -> W -> ...` and gives `X -> W` its part of the proportional split of its own length. Each later child attaches at the end of that chain and keeps only the tail of its split (`fn reroute()` [packages/treetime-graph/src/pass_through.rs#L55-L92](../../packages/treetime-graph/src/pass_through.rs#L55-L92), described in [kb/algo/graph.md](../algo/graph.md#merging-trees-by-name)). The merge visits the edges of a segment in the order of the child names (`SegmentTree::parent` is a `BTreeMap`), so B's length of `X -> W` depends on node names. The graph then gives the shared edge the mean over its segments ([kb/decisions/segment-graph-shared-edge-mean-length.md](../decisions/segment-graph-shared-edge-mean-length.md)).

Only `ancestral` with the input branch lengths uses these lengths. `optimize` and `timetree` fit the lengths again and use them only as the start.

> [!IMPORTANT]
> **Decision required.** Options for B's length of `X -> W`:
>
> - The part of the first child in name order (current)
> - The mean of the head parts of all children of `X` that move below `W`

## Validation

A merge test with a polytomy whose children give different head parts: the length of `X -> W` does not change when the children are renamed.
