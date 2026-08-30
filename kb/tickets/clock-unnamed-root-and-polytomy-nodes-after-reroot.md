# Name the rerooted root and polytomy nodes in the clock command

After rerooting, the new root node (and any nodes created during polytomy resolution) have no name. v0 assigns `NODE_XXXXXXX` names to unnamed internal nodes after every reroot via `_prepare_nodes()` (`treeanc.py:471-478`), called from `prepare_tree()` at the end of the `reroot()` method. v1's `assign_node_names` (`packages/treetime-graph/src/assign_node_names.rs:7`) runs during Newick parsing (`nwk.rs:100`) but is not called by the reroot path itself.

The `timetree` command already names every unnamed node once after its pipeline (`packages/treetime/src/commands/timetree/run.rs`). Apply the equivalent to `clock`.

## Current state

`create_new_root_node` (`packages/treetime/src/clock/reroot.rs:163`) creates a new node via `N::default()` which has `name=None`. In the `clock` command this node stays unnamed through output. The auspice writer assigns a fallback name `node_<key>` for display, but other consumers (TSV output, Newick annotations) see the unnamed state.

## Task

Call `assign_node_names` in the `clock` command after rerooting completes, or move the call into the shared reroot path so every rerooting command names its new nodes uniformly, matching v0's `prepare_tree()` pattern. Prefer the shared path if it does not regress the `timetree` command's existing post-pipeline call.

## Impact

- Clock `confidence_intervals.tsv` includes unnamed nodes with an empty name column
- Clock Newick/Nexus output uses default formatting for the unnamed root
- Differs from v0, where all internal nodes have `NODE_` names after any topology change

## Related issues

- Source: [kb/issues/N-clock-unnamed-root-after-reroot.md](../issues/N-clock-unnamed-root-after-reroot.md) -- delete after full resolution
