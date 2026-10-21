# Unnamed nodes are named differently by each writer

Input parsing names unnamed nodes with `assign_node_names()` ([packages/treetime-graph/src/assign_node_names.rs#L7](../../packages/treetime-graph/src/assign_node_names.rs#L7)), and several pipelines call it again later (timetree, optimize, prune). Nodes created during a run, such as the new root of a `clock` reroot, can stay unnamed, and each writer then invents its own name: an empty string in UShER MAT ([packages/app-output/src/usher_mat.rs#L32](../../packages/app-output/src/usher_mat.rs#L32)), `node_<key>` in Auspice JSON, augur node data, the mugration CSVs and the amino-acid FASTA (`node_name_or_key()`, [packages/treetime-graph/src/assign_node_names.rs#L66-L68](../../packages/treetime-graph/src/assign_node_names.rs#L66-L68)), no label in Newick, and `NODE_0000001` in the reconstructed FASTA. The same node gets different names in the files of one run.

Required behavior: the project convention is `NODE_0000001` (7 digits, numbered in preorder, skipping names already in use), as v0 names nodes in `_prepare_nodes()` ([packages/legacy/treetime/treetime/treeanc.py#L460-L479](../../packages/legacy/treetime/treetime/treeanc.py#L460-L479)). Naming happens outside the formats: during input parsing for nodes that exist in the input, and at node creation for nodes created during a run (new root, polytomy resolution, merged nodes). Every node then has a name, the name map holds `String` instead of `Option<String>`, the repeated `assign_node_names()` calls go away, and writers never invent names.

## Generated names differ from v0

`assign_node_names()` gives a node the same generated name as v0 only on a ladderized tree whose internal nodes are all unnamed. Two differences of the numbering ([packages/treetime-graph/src/assign_node_names.rs#L25-L39](../../packages/treetime-graph/src/assign_node_names.rs#L25-L39) against v0 `prepare_tree()` and `_prepare_nodes()`, [packages/legacy/treetime/treetime/treeanc.py#L446-L480](../../packages/legacy/treetime/treetime/treeanc.py#L446-L480)):

- **Node order**: v0 ladderizes the tree before it names the nodes (unless `keep_node_order` is set), and v1 names them in the preorder of the input file. On a tree that is not ladderized, the same `NODE_*` name points to different nodes in v0 and v1 outputs
- **Counter**: v0 advances the counter at every internal node, also at an internal node that has a name in the input; v1 advances it only when it assigns a name. On a tree with named internal nodes, the numbers that follow differ

The mugration golden-master test names the v1 nodes as v0 does to compare them ([packages/treetime/src/mugration/__tests__/test_gm_mugration.rs](../../packages/treetime/src/mugration/__tests__/test_gm_mugration.rs)); the outputs of every command keep the v1 names. On `dengue`, `tb`, `rsv`, and `mpox` (20 tips), the same `NODE_*` name in v0 and v1 names different nodes.

> [!IMPORTANT]
> **Decision required.** v0 parity of generated names conflicts with the output order of the tree: v1 writes the tree in input order unless a tree-ordering flag such as `--ladderize` is given (`TopologyOrderArgs`, [packages/app-commands/src/commands/shared/topology_order_args.rs](../../packages/app-commands/src/commands/shared/topology_order_args.rs)), and v0 always ladderizes. With v0 numbering, the names of a tree written in input order do not follow its preorder. Options:
>
> - Number the nodes as v0 does (ladderized preorder, counter advanced at every internal node), independent of the output order
> - Number the nodes in the preorder of the tree as written, and record the difference from v0

## Related issues

- [N-clock-unnamed-root-after-reroot.md](N-clock-unnamed-root-after-reroot.md)
- [M-io-confidence-tsv-unnamed-root-has-empty-name.md](M-io-confidence-tsv-unnamed-root-has-empty-name.md)
