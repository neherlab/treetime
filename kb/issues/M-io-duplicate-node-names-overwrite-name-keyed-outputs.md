# Duplicate node names overwrite entries in name-keyed outputs

An input tree may name two nodes the same, for example two leaves `A`. Every run keeps both names, and the outputs that key their entries by node name then keep one entry per name without a message. The data of the other node is lost from these files, while the tree files still contain both nodes.

## Behavior

- `fn assign_node_names()` [packages/treetime-graph/src/assign_node_names.rs#L7](../../packages/treetime-graph/src/assign_node_names.rs#L7) keeps every input name, duplicates included, and names only the unnamed nodes (`NODE_0000000`, ...)
- Augur node data writes `nodes` as a map from name to fields. A later node replaces an earlier node of the same name in all three layouts: `refine` (timetree, optimize) [packages/app-output/src/augur_node_data_refine.rs#L27-L55](../../packages/app-output/src/augur_node_data_refine.rs#L27-L55), `ancestral` [packages/app-output/src/augur_node_data_ancestral.rs#L67](../../packages/app-output/src/augur_node_data_ancestral.rs#L67), and `traits` (mugration) [packages/app-output/src/augur_node_data_traits.rs#L90](../../packages/app-output/src/augur_node_data_traits.rs#L90). The augur format is a JSON object keyed by name, so it cannot hold two nodes of the same name
- The mugration traits CSV collects its rows into an `IndexMap` keyed by name [packages/app-output/src/trait_tables.rs#L28-L36](../../packages/app-output/src/trait_tables.rs#L28-L36): it writes one row per name, at the position of the first node and with the value of the last
- The mugration confidence CSV writes one row per node, so duplicate names give two rows with the same name [packages/app-output/src/trait_tables.rs#L48-L56](../../packages/app-output/src/trait_tables.rs#L48-L56)
- The input check lists duplicate leaf names in `duplicate_tip_names` [packages/app-commands/src/check_inputs.rs#L328](../../packages/app-commands/src/check_inputs.rs#L328), but no analysis command warns about them

## v0 behavior

v0 has no duplicate check either. `TreeAnc` looks leaves up in a dictionary keyed by name, which keeps the last leaf of each name (`_leaves_lookup`, [packages/legacy/treetime/treetime/treeanc.py#L459](../../packages/legacy/treetime/treetime/treeanc.py#L459)), and the mugration `confidence.csv` writes one row per node ([packages/legacy/treetime/treetime/wrappers.py#L909-L910](../../packages/legacy/treetime/treetime/wrappers.py#L909-L910)). v0 writes no augur node data.

## Required behavior

[kb/decisions/duplicate-names-warned-ids-from-input-order.md](../decisions/duplicate-names-warned-ids-from-input-order.md) sets the rule: warn on every surface, keep running, identify nodes by input order, and never guess which node a name means.

- Warn about duplicate node names on every surface. The input check finds them in linear time, but no command reports them
- Write the traits CSV with one row per node in node order, as the confidence CSV does, because a CSV file can hold repeated names
- Keep writing augur node data: its `nodes` object cannot hold two nodes of the same name, so the warning is the only remedy

## Validation

- A tree with two leaves `A` gives the warning on every surface, two rows `A` in both mugration CSVs, and one entry `A` in augur node data of every layout
- A tree with unique names produces byte-identical outputs
