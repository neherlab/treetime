# Rerooted root and polytomy-resolution nodes stay unnamed in the clock command

> [!WARNING]
> **Needs review.** The output impact list below is not re-verified against the current writers. Tree writers fall back to `node_<key>` for unnamed nodes (`packages/app-output/src/tree_output.rs:659-661`), so Newick and Nexus output may already show a fallback name rather than an empty one. The clock table takes the raw name (`packages/treetime/src/clock/rtt.rs:26`), so its row for the new root has no name.

After rerooting, the new root node (and any nodes created during polytomy resolution) have no name. v0 assigns `NODE_XXXXXXX` names to unnamed internal nodes after every reroot via `_prepare_nodes()` (`treeanc.py:471-478`), called from `prepare_tree()` at the end of the `reroot()` method. v1's `assign_node_names` (`packages/treetime-graph/src/assign_node_names.rs:7`) runs during Newick parsing (`packages/treetime-io/src/nwk.rs:84`) but is not called by the reroot path itself, so a command that reroots without a later naming pass emits the unnamed root.

The `timetree` command does not have this defect: it names every unnamed node after the pipeline, before serialization (`packages/app-commands/src/commands/timetree/run.rs:362`). The `clock` command reroots (`packages/app-commands/src/commands/clock/run.rs`) and has no post-reroot naming pass.

## Current state

`create_new_root_node` (`packages/treetime/src/clock/reroot.rs:133`) splits the root edge with `split_edge`, which adds a node that has no entry in the name map. The clock pipeline then calls `restrict_node_names` (`packages/treetime/src/clock/pipeline.rs:67`), which keeps the names of existing nodes and gives new nodes `None`. In the `clock` command this node stays unnamed through output.

## Impact

- The clock table (`ClockCsv`) has an empty name for the new root
- Tree outputs show at most the generic `node_<key>` fallback for the new root
- Differs from v0, where all internal nodes have `NODE_` names after any topology change

## Fix

Give the `clock` command a post-reroot naming pass, or move `assign_node_names` into the shared reroot path so every rerooting command names its new nodes, matching v0's `prepare_tree()` pattern. Prefer the shared path if it does not regress the `timetree` command's existing post-pipeline naming.
