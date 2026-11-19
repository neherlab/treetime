# A shared edge of a segment graph starts from the mean of its input lengths

When several segment trees contain the same parent-child pair, the merged graph has one edge with one branch length: the mean of the lengths that these trees give. A tree without a length for the edge does not count, and an edge that no tree gives a length has none. Segment trees estimated separately give a shared edge different lengths, and the merge does not check them against each other ([kb/algo/graph.md](../algo/graph.md#merging-trees-by-name)).

**Type**: Scientific default of a new v1 input (the starting branch lengths of a segment graph). Approved.

Code: `fn merge_trees_by_name` [packages/treetime-graph/src/merge_by_name.rs](../../packages/treetime-graph/src/merge_by_name.rs).

## Effect on results

- **`ancestral` with input branch lengths**: every segment is reconstructed with the mean length on its shared edges, not with the length of its own tree. Where the lengths of the trees differ, a segment's sequences, mutations and log-likelihood can differ from a run on its own tree. The toy graph of the tests has trees that agree on every shared length, so its segments match their own runs exactly ([packages/treetime/src/ancestral/__tests__/test_pipeline_segment_graph.rs](../../packages/treetime/src/ancestral/__tests__/test_pipeline_segment_graph.rs))
- **`optimize` and `timetree`**: they fit one length per shared edge to all segments that carry it, so the mean is only the starting point ([kb/reports/2026-11-15-segment-validation-optimize.md](../reports/2026-11-15-segment-validation-optimize.md))
- **Pass-through nodes**: when the merge routes a segment's edge through nodes of another segment, it splits that segment's length over the new edges in proportion to the other segment's lengths before the mean is taken ([kb/algo/graph.md](../algo/graph.md#merging-trees-by-name))

## Why the mean

- **No segment is privileged**: the lengths of the trees are separate estimates of the same branch, and the mean pools them
- **Shared lengths are the model**: an edge has one length for all segments it carries, as one branch of a tree has one length for all sites ([kb/proposals/edge-scoped-inheritance.md](../proposals/edge-scoped-inheritance.md))
