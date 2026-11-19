# The merge of segment trees moves an edge past a reassortment when the trees hold every pass-through node

`fn insert_pass_through_nodes()` [packages/treetime-graph/src/pass_through.rs#L102-L154](../../packages/treetime-graph/src/pass_through.rs#L102-L154) moves a segment's edge `X -> Y` below the nodes `W1 ... Wk` that another segment puts on its path from `X` to `Y` ([kb/algo/graph.md](../algo/graph.md#merging-trees-by-name)). This assumes that the segment's lineage between `X` and `Y` is the lineage of the other segment. Trees from tree inference omit the nodes that a segment passes through, so a direct edge in such a tree usually means that the segment passes through `W` without branching there, and the assumption holds unless a reassortment on that stretch went undetected.

A tree that holds every node that its segment passes through (the true local trees of a simulation, or an ARG) states with a direct edge `X -> Y` that the segment bypasses `W`: a reassortment separates the segment from `W`. The merge cannot tell the two kinds of input apart, so it moves the edge below `W` and deletes a true parent edge of `Y`.

## Evidence

[kb/reports/2026-11-15-segment-validation-timetree.md](../reports/2026-11-15-segment-validation-timetree.md) (section "Merge of trees that hold all their pass-through nodes"): on the true local trees of 6 simulated replicates, the merge reroutes one edge in the `full` form (`node_0132 -> node_0126` of `high/rep00`) and one in the `kept` form (`node_0017 -> node_0003` of `high/rep02`). The same step makes contraction of hidden nodes in the merged graph differ from contraction of the input trees by one edge in 1 of 6 replicates (`high/rep00`, section "Contraction of the input and contraction in the graph"). Collapsed trees, the input that tree inference gives, are not affected. The effect of the rerouted edges on the dates was not measured apart from the other differences between the forms.

> [!IMPORTANT]
> **Decision required.** The merge needs information that topology does not give. Options:
>
> - The input states that its trees hold every pass-through node (a flag, or an ARG file as input), and the merge then inserts no nodes
> - TreeKnit gives the attachment points of reassortant lineages, which give the pass-through nodes directly ([M-io-treeknit-output-reaches-segment-graphs-only-through-a-python-prototype.md](M-io-treeknit-output-reaches-segment-graphs-only-through-a-python-prototype.md))
> - Keep the assumption as a documented limitation of the merge

## Validation

`validation/segments.py evaluate` lists the rerouted edges in the table "Merged graphs" of `summary.md`; with the fix, the `full` and `kept` forms have no rerouted edges, and the `collapsed` form keeps its dates.
