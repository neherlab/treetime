# The new root after rerooting stays unnamed in the clock command

After rerooting, the new root node has no name. v0 assigns `NODE_XXXXXXX` names to unnamed internal nodes after every reroot via `_prepare_nodes()` (`treeanc.py:471-478`), called from `prepare_tree()` at the end of the `reroot()` method. v1's `fn assign_node_names()` [packages/treetime-graph/src/assign_node_names.rs#L7](../../packages/treetime-graph/src/assign_node_names.rs#L7) runs during Newick parsing [packages/treetime-io/src/nwk.rs#L90](../../packages/treetime-io/src/nwk.rs#L90) but is not called by the reroot path itself, so a command that reroots without a later naming pass emits the unnamed root.

The `timetree` command does not have this defect: it names every unnamed node after polytomy resolution [packages/treetime/src/timetree/round.rs#L169](../../packages/treetime/src/timetree/round.rs#L169) and before the final outputs [packages/treetime/src/timetree/pipeline.rs#L474](../../packages/treetime/src/timetree/pipeline.rs#L474). The `clock` command reroots and has no naming pass after the reroot.

## Current state

`fn reroot_at_best_root()` [packages/treetime/src/clock/reroot.rs#L41](../../packages/treetime/src/clock/reroot.rs#L41) applies the new root through `fn apply_root_at_edge()` [packages/treetime/src/reroot/orchestrate.rs#L19](../../packages/treetime/src/reroot/orchestrate.rs#L19), which splits the chosen edge with `split_edge` of `treetime_graph::reroot` [packages/treetime/src/reroot/orchestrate.rs#L50](../../packages/treetime/src/reroot/orchestrate.rs#L50). The new node has no entry in the name map. The clock pipeline then calls `fn restrict_node_names()` [packages/treetime/src/clock/pipeline.rs#L56](../../packages/treetime/src/clock/pipeline.rs#L56), which keeps the names of existing nodes and gives new nodes `None`. In the `clock` command this node stays unnamed through output.

## Impact

`treetime clock --tree=data/flu/h3n2/20/tree.nwk --metadata=data/flu/h3n2/20/metadata.tsv --output-all=<dir>` with default flags gives:

- `clock.clock.csv`: the row of the new root has an empty name (`,0.0,,1996.75...,,`), because `fn gather_clock_regression_results()` [packages/treetime/src/clock/rtt.rs#L27](../../packages/treetime/src/clock/rtt.rs#L27) copies the raw name
- `clock.nwk` and `clock.nexus`: the root has no label (the tree ends in `...:0.00408);`), while every other internal node has a `NODE_` name

This differs from v0, where all internal nodes have `NODE_` names after any topology change.

## Fix

Give the `clock` command a post-reroot naming pass, or move `assign_node_names` into the shared reroot path so every rerooting command names its new nodes, matching v0's `prepare_tree()` pattern. Prefer the shared path if it does not regress the existing naming of the `timetree` command.
