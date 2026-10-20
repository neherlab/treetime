# Unnamed nodes are named differently by each writer

Input parsing names unnamed nodes with `assign_node_names()` ([packages/treetime-graph/src/assign_node_names.rs#L7](../../packages/treetime-graph/src/assign_node_names.rs#L7)), and several pipelines call it again later (timetree, optimize, prune). Nodes created during a run, such as the new root of a `clock` reroot, can stay unnamed, and each writer then invents its own name: an empty string in UShER MAT ([packages/app-output/src/usher_mat.rs#L32](../../packages/app-output/src/usher_mat.rs#L32)), `node_<key>` in Auspice JSON, augur node data, the mugration CSVs and the amino-acid FASTA (`node_name_or_key()`, [packages/treetime-graph/src/assign_node_names.rs#L66-L68](../../packages/treetime-graph/src/assign_node_names.rs#L66-L68)), no label in Newick, and `NODE_0000001` in the reconstructed FASTA. The same node gets different names in the files of one run.

Required behavior: the project convention is `NODE_0000001` (7 digits, numbered in preorder, skipping names already in use), as v0 names nodes in `_prepare_nodes()` ([packages/legacy/treetime/treetime/treeanc.py#L460-L479](../../packages/legacy/treetime/treetime/treeanc.py#L460-L479)). Naming happens outside the formats: during input parsing for nodes that exist in the input, and at node creation for nodes created during a run (new root, polytomy resolution, merged nodes). Every node then has a name, the name map holds `String` instead of `Option<String>`, the repeated `assign_node_names()` calls go away, and writers never invent names.

## Related issues

- [N-clock-unnamed-root-after-reroot.md](N-clock-unnamed-root-after-reroot.md)
- [M-io-confidence-tsv-unnamed-root-has-empty-name.md](M-io-confidence-tsv-unnamed-root-has-empty-name.md)
