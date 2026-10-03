# Input branch support stays on its node when a reroot inverts the branch

A Newick support value describes the split of the branch above a node. v1 keys the parsed values by node (`fn NwkParse::confidences()` in [packages/treetime-io/src/nwk.rs](../../packages/treetime-io/src/nwk.rs)), and a reroot inverts the edges on the path between the old and the new root and keeps every node key (`fn apply_reroot_topology()` in [packages/treetime-graph/src/reroot.rs](../../packages/treetime-graph/src/reroot.rs)). No step moves the support values, so after a reroot every node on that path reports the support of a split it no longer sits below.

## Affected outputs

- `treetime timetree` augur node data: the per-node `confidence` field (`fn write_augur_node_data_json()` in [packages/app-output/src/augur_node_data.rs](../../packages/app-output/src/augur_node_data.rs)). Timetree reroots by default
- `treetime optimize --reroot`: the Auspice node `confidence` and the augur node data `confidence`
- Outputs of runs that keep the root are correct

`timetree` Auspice JSON does not write input support values until this is fixed, while `ancestral`, `optimize`, `prune` and `mugration` Auspice JSON do.

## Fix direction

Key input support by edge, or move each value to the inverted edge's new child when a reroot inverts an edge, so the value stays with its split. Then write it in the timetree Auspice JSON as the other commands do.
