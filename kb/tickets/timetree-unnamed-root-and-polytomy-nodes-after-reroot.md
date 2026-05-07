# Rerooted root and polytomy-resolution nodes stay unnamed

After rerooting, the new root node (and any nodes created during polytomy resolution) have no name. v0 assigns `NODE_XXXXXXX` names to unnamed internal nodes after every reroot via `_prepare_nodes()` (`treeanc.py:471-478`), called from `prepare_tree()` at the end of the `reroot()` method. v1's equivalent `assign_node_names` (`treetime-graph/src/assign_node_names.rs:7`) only runs during Newick parsing (`nwk.rs:100`) and is never called after rerooting or topology changes.

## Current state

`create_new_root_node` (`packages/treetime/src/commands/clock/reroot.rs:163`) creates a new node via `N::default()` which has `name=None`. This node stays unnamed throughout the pipeline. The auspice writer assigns a fallback name `node_<key>` for display, but other consumers (TSV output, Newick annotations) see the unnamed state.

## Impact

- `confidence_intervals.tsv` includes unnamed nodes with empty name column
- Newick/Nexus output uses default formatting for unnamed nodes
- Differs from v0 where all internal nodes have `NODE_` names after any topology change

## Fix

Call `assign_node_names` after rerooting and polytomy resolution, matching v0's `prepare_tree()` pattern. This is a broader fix that affects all commands using reroot (timetree, clock).

## Related issues

- Source: [N-timetree-unnamed-root-after-reroot.md](../issues/N-timetree-unnamed-root-after-reroot.md) -- delete after full resolution
