# Unnamed internal nodes report a stale root-to-node divergence

The timetree output column `div` (`TimetreeNodeOut.div`) is computed once, after the last time inference, by `fn final_divergences` ([packages/treetime/src/timetree/divergence.rs#L9](../../packages/treetime/src/timetree/divergence.rs#L9)), called from `fn gather_results` of the timetree pipeline ([packages/treetime/src/timetree/pipeline.rs#L444](../../packages/treetime/src/timetree/pipeline.rs#L444)). The value of a node depends on whether it has a name in the node names of the last time inference (the names before the final `assign_node_names`):

- Named nodes get the root-to-node divergence on the final tree and final branch lengths
- Unnamed nodes get the divergence the clock filter computed for them in `fn filter_clock_outliers` ([packages/treetime/src/timetree/pre_loop.rs#L243](../../packages/treetime/src/timetree/pre_loop.rs#L243)), when the filter ran and the node existed at that time
- All other unnamed nodes get `0.0`: every unnamed node when the clock filter did not run (`--clock-filter=0`), and nodes that a reroot created later (the split node of the post-ancestral reroot)

The clock-filter value is computed on the topology and branch lengths at the time of the filter, which is before the ML post-step (marginal mode) and before the post-ancestral reroot. Both change branch lengths, and the reroot changes the root, so the reported value of an unnamed internal node does not describe the final tree.

Polytomy resolution names every node when it changes the tree, so after such a round all nodes take the first rule.

## Impact

The `div` attribute of unnamed internal nodes in the timetree tree outputs (Newick comments, Nexus, Auspice JSON, augur node data) is stale or `0.0`, while named nodes in the same tree carry current values. Input trees whose internal nodes carry names do not show the effect. Named leaves always carry current values.

## Cause

`final_divergences` refreshes the divergence only for nodes that have a name in the node names of the last time inference, and falls back to the clock-filter value for the others.

## Open question

Decide whether every node of the final tree reports the root-to-node divergence on the final tree, independent of its name.
