# Clock Newick output fuses input branch support into node names

## v0 location

`estimate_clock_model()` writes the tree of `treetime clock` with `Phylo.write(myTree.tree, outtree_name, 'newick')` and does not clear the input support values first [packages/legacy/treetime/treetime/wrappers.py#L1032-L1041](../../packages/legacy/treetime/treetime/wrappers.py#L1032-L1041). The file is `rerooted.newick` by default. With `--keep-root`, it is `pruned.newick` when `--prune-outliers` is also given and `.output.newick` otherwise.

## Erratum

A Newick support value (bootstrap, SH-like local support, posterior probability) belongs to the split of the branch above its node. Tree builders write it as the label of an internal node, for example `(C,(D,E)70)80`. Biopython's reader moves a label that parses as a number into `clade.confidence` and clears the node name [[src](https://github.com/biopython/biopython/blob/d7e4b8b19399668b09442a5b35765d9186b5f665/Bio/Phylo/NewickIO.py#L231-L241)]. `TreeAnc._prepare_nodes()` then gives every unnamed internal node a name, `NODE_` and seven digits [packages/legacy/treetime/treetime/treeanc.py#L469-L478](../../packages/legacy/treetime/treetime/treeanc.py#L469-L478). The clock output has two defects.

- **Fused labels**: Biopython's Newick writer writes the node name and then the confidence with the format `%1.2f`, with no separator [[src](https://github.com/biopython/biopython/blob/d7e4b8b19399668b09442a5b35765d9186b5f665/Bio/Phylo/NewickIO.py#L296-L308)] [[src](https://github.com/biopython/biopython/blob/d7e4b8b19399668b09442a5b35765d9186b5f665/Bio/Phylo/NewickIO.py#L364-L378)]. The node `NODE_0000016` with support `1.0` is written as `NODE_00000161.00:0.017777875`. This label is neither the node name nor a support value. A reader takes the whole label as a name, so the support values are lost, and the third decimal of each value is cut off (`0.959` is written as `0.96`)
- **Values on the wrong split**: Biopython's `Tree.root_with_outgroup()` moves the branch lengths along the path between the old and the new root, but each node keeps its `confidence` [[src](https://github.com/biopython/biopython/blob/d7e4b8b19399668b09442a5b35765d9186b5f665/Bio/Phylo/BaseTree.py#L873-L878)]. After a reroot, every node on that path carries the value of another split. The clock filter reroots also when `--keep-root` is given ([clock-keep-root-ignored-by-clock-filter.md](clock-keep-root-ignored-by-clock-filter.md))

## Evidence

- Adjacent code clears the values: the shared export of `ancestral`, `timetree` and `arg` sets `n.confidence = None` before it writes Newick and Nexus [packages/legacy/treetime/treetime/CLI_io.py#L160](../../packages/legacy/treetime/treetime/CLI_io.py#L160), and `mugration` does the same [packages/legacy/treetime/treetime/wrappers.py#L897](../../packages/legacy/treetime/treetime/wrappers.py#L897). The clock writer is the only v0 tree writer that keeps them
- The output contradicts the label rule of v0's own reader: Biopython reads a label as support only when the whole label is a number [[src](https://github.com/biopython/biopython/blob/d7e4b8b19399668b09442a5b35765d9186b5f665/Bio/Phylo/NewickIO.py#L66-L76)]. Biopython 1.88 reads no support values from the clock output of `data/flu/h3n2/20`
- `treetime clock --tree data/flu/h3n2/20/tree.nwk --dates data/flu/h3n2/20/metadata.tsv --sequence-length 1400` writes 16 internal labels with a fused value. 9 of them sit on a split whose input value is different. With `--keep-root --clock-filter 0`, the root does not move and all 16 values sit on their own split

## v0 impact

- The support values of the input tree cannot be read back from the clock output
- The internal node names of the clock tree differ from the names in `rtt.csv` of the same run (`NODE_00000161.00` against `NODE_0000016`), so the two files cannot be joined on internal nodes

## v1 status

v1 `clock` writes plain node names and no support values in its Newick and Nexus output; nodes without a name get `NODE_` names from `assign_node_names()` [packages/treetime-graph/src/assign_node_names.rs](../../packages/treetime-graph/src/assign_node_names.rs). On `data/flu/h3n2/20`, `clock.nwk` has labels such as `NODE_0000013`. No v1 output writes input support values; their placement after a reroot and their output format are tracked in [kb/issues/M-io-branch-support-dropped-from-outputs.md](../issues/M-io-branch-support-dropped-from-outputs.md).
