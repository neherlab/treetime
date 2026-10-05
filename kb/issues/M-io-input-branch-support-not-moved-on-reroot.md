# Input branch support stays on its node when a reroot inverts the branch

A Newick support value describes the split of the branch above a node. v1 keys the parsed values by node (`fn NwkParse::confidences()` in [packages/treetime-io/src/nwk.rs](../../packages/treetime-io/src/nwk.rs)), and a reroot inverts the edges on the path between the old and the new root and keeps every node key (`fn apply_reroot_topology()` in [packages/treetime-graph/src/reroot.rs](../../packages/treetime-graph/src/reroot.rs)). No step moves the support values, so after a reroot every node on that path reports the support of a split it no longer sits below.

Example: the input root has children `A`, `B` and `X = (C, (D, E)70)80`. After a reroot on the branch to `E`, the node above `{A, B, C}` shows 80, although 80 was measured for `{A, B}` against `{C, D, E}`; its correct value is 70, and the node above `{A, B}` shows nothing instead of 80.

## Affected outputs

- `treetime timetree` augur node data: the per-node `confidence` field (`fn write_augur_node_data_json()` in [packages/app-output/src/augur_node_data.rs](../../packages/app-output/src/augur_node_data.rs)). Timetree reroots by default
- `treetime optimize --reroot`: the Auspice node `confidence` and the augur node data `confidence`
- `treetime clock` Auspice JSON: writes no input support values (`fn clock_to_auspice()` in [packages/app-output/src/clock_tree_output.rs](../../packages/app-output/src/clock_tree_output.rs)), and `ClockNodeOut` has no field for them. Clock reroots by default, so adding the values needs the same remap
- Outputs of runs that keep the root are correct

`timetree` and `clock` Auspice JSON do not write input support values until this is fixed, while `ancestral`, `optimize`, `prune` and `mugration` Auspice JSON do.

## Decided placement rules

Support is stored per edge: the Newick parse keys each value by the edge above its node and drops a value on the root, which has no branch. gotree stores support on the edge object, so a reroot that only flips edge direction keeps each value with its split (`Reroot()` and `ReorderEdges()` in `tree/tree.go` of gotree). Each topology edit then follows one rule (approved 2026-10-05):

- An inverted edge keeps its key, so its value stays
- A new root placed inside a branch (`split_edge()` in [packages/treetime-graph/src/reroot.rs](../../packages/treetime-graph/src/reroot.rs)): both halves get the branch's value
- Removing a node with one parent and one child (`remove_node_if_trivial()`: the old root after a reroot, polytomy cleanup): the merged branch gets the larger of the two values, or none when it ends at a leaf, as gotree does when it removes a two-child root (`UnRoot()` in `tree/tree.go`)
- Removing a one-child root (`remove_stem_root()`): the removed branch's value is dropped
- Collapsing a branch (`collapse_edge()` in [packages/treetime/src/optimize/topology/collapse.rs](../../packages/treetime/src/optimize/topology/collapse.rs)): its value is dropped; the child branches keep their keys and values
- A branch created for a new split (polytomy resolution in [packages/treetime/src/timetree/optimization/polytomy/apply.rs](../../packages/treetime/src/timetree/optimization/polytomy/apply.rs)): no value
- Reordering children (`topology_order.apply()`) keeps edge keys, so values stay

Once support follows its split, `timetree` and `clock` write it in Auspice JSON as the other commands do. v0 has no counterpart: it writes no input support after a timetree run (`n.confidence = None`, [packages/legacy/treetime/treetime/CLI_io.py#L160](../../packages/legacy/treetime/treetime/CLI_io.py#L160)), and its Auspice `confidence` is a substitute from mutation counts, `1 - exp(-n)` ([packages/legacy/treetime/treetime/CLI_io.py#L321-L327](../../packages/legacy/treetime/treetime/CLI_io.py#L321-L327)). Only v0 `clock` writes input support, fused into the node names of its Newick output and on the wrong splits after a reroot ([kb/v0-errata/clock-newick-support-fused-into-node-names.md](../v0-errata/clock-newick-support-fused-into-node-names.md)).

## Open: how support travels through a run

The pipelines that edit topology (`clock`, `timetree`, `optimize`, `prune`) do not see support today; the runners hold it keyed by node. Timetree, for example, reroots and then resolves polytomies, and each step adds, removes or merges edges. Options:

- Pipelines take the per-edge support map as input, update it where they update branch lengths for a topology edit (`record_split()`, `record_merge()` and the other edit sites), and return it. Small and follows the branch-length pattern; the core carries a value it never computes with
- Pipelines return a log of their topology edits, and the runners replay it on the support map. The core never sees support; every edit, including the polytomy merges inside timetree's iterations, must be recorded exactly
- One per-edge record holds branch length and support, and every edit updates it once. The cleanest model; it changes every branch-length user in the core
