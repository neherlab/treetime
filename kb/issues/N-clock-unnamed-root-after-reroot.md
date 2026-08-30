# Rerooted root and polytomy-resolution nodes stay unnamed in the clock command

After rerooting, the new root node (and any nodes created during polytomy resolution) have no name. v0 assigns `NODE_XXXXXXX` names to unnamed internal nodes after every reroot via `_prepare_nodes()` (`treeanc.py:471-478`), called from `prepare_tree()` at the end of the `reroot()` method. v1's `assign_node_names` (`packages/treetime-graph/src/assign_node_names.rs:7`) runs during Newick parsing (`nwk.rs:100`) but is not called by the reroot path itself, so a command that reroots without a later naming pass emits the unnamed root.

The `timetree` command no longer has this defect: it names every unnamed node once after the pipeline, before serialization (`packages/treetime/src/commands/timetree/run.rs`, immediately after `pipeline::run`). The `clock` command reroots (`packages/treetime/src/commands/clock/run.rs`) and still has no post-reroot naming pass.

## Current state

`create_new_root_node` (`packages/treetime/src/clock/reroot.rs:163`) creates a new node via `N::default()` which has `name=None`. In the `clock` command this node stays unnamed through output. The auspice writer assigns a fallback name `node_<key>` for display, but other consumers (TSV output, Newick annotations) see the unnamed state.

## Impact

- Clock `confidence_intervals.tsv` includes unnamed nodes with an empty name column
- Clock Newick/Nexus output uses default formatting for the unnamed root
- Differs from v0, where all internal nodes have `NODE_` names after any topology change

## Fix

Give the `clock` command a post-reroot naming pass, or move `assign_node_names` into the shared reroot path so every rerooting command names its new nodes, matching v0's `prepare_tree()` pattern.

## Related tickets

- [kb/tickets/clock-unnamed-root-and-polytomy-nodes-after-reroot.md](../tickets/clock-unnamed-root-and-polytomy-nodes-after-reroot.md)
