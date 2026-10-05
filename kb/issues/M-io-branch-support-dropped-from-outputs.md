# Branch support values are dropped from every output

v1 writes no branch support in any output. It reads support values from the input Newick and discards them, and it computes no support of its own. This covers two kinds of value:

- **Input support**: values that a tree builder wrote into the input tree, for example bootstrap percentages
- **Mutation-based support**: the substitute that v0 computes from the number of mutations on a branch and writes into its Auspice JSON

Both are left out on purpose until the open questions below are decided (approved 2026-10-05). Input support needs a placement rule that survives the topology edits of a run. The meaning of the Auspice `confidence` attribute needs a decision, because v0 and augur write different quantities under that key.

## Background

### What a support value is

A support value measures how strongly the data support one **split**: the division of the leaves that a branch makes, for example `{A, B}` against `{C, D, E}`. A split does not depend on the root, so a support value belongs to a branch of the unrooted tree, not to a node. Common kinds, on different scales:

- **Bootstrap** ([Felsenstein 1985](https://doi.org/10.2307/2408678) [[1](#ref-1)]): the fraction of trees, rebuilt from resampled alignment columns, that contain the split. 0-100 or 0-1
- **UFBoot** ([Minh et al. 2013](https://doi.org/10.1093/molbev/mst024) [[2](#ref-2)]): a fast bootstrap approximation of IQ-TREE (`-B`). 0-100, and not comparable to the standard bootstrap
- **SH-aLRT and SH-like local support** ([Anisimova & Gascuel 2006](https://doi.org/10.1080/10635150600755453) [[3](#ref-3)]): a local likelihood test of the branch against the two other arrangements around it. FastTree computes SH-like local support by default ([Price et al. 2010](https://doi.org/10.1371/journal.pone.0009490) [[4](#ref-4)]) and writes it as 0-1 with three decimals
- **Posterior probability** (MrBayes, BEAST): usually in comments such as `[&posterior=0.98]`, not in labels
- **Transfer bootstrap expectation** ([Lemoine et al. 2018](https://doi.org/10.1038/s41586-018-0043-0) [[5](#ref-5)]): a bootstrap variant for large trees

A support value alone does not tell which method produced it. Readers use it to judge which clades to trust, viewers color or label branches by it, and tools collapse weakly supported branches into polytomies before further analysis (`collapse support` of gotree).

### How Newick stores it

Newick has no field for branch data. Tree builders write the value as the label of the node below the branch: `(C,(D,E)70)80`. A reader must guess: Biopython ([[src](https://github.com/biopython/biopython/blob/d7e4b8b19399668b09442a5b35765d9186b5f665/Bio/Phylo/NewickIO.py#L231-L241)]) and v1 (`fn visit_internal()` [packages/util-newick/src/parse.rs#L151-L157](../../packages/util-newick/src/parse.rs#L151-L157)) both read an internal label as support when the whole label parses as a number, and as a name otherwise.

Because the format stores the value on a node, tools that reroot a tree often leave the value on its node and so on the wrong split. [Czech et al. 2017](https://doi.org/10.1093/molbev/msx055) [[6](#ref-6)] found this defect in many tree viewers and toolkits.

The label rule has two traps:

- IQ-TREE run with several support methods writes them in one label separated by `/`, for example `80.5/95` ([[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/main/phyloanalysis.cpp#L1100-L1115)]). The label is not a number, so v1 keeps it as the node name
- `fn parse_support_value()` uses `f64::from_str` ([packages/util-newick/src/parse.rs#L275-L277](../../packages/util-newick/src/parse.rs#L275-L277)), which accepts `inf` and `nan`. An internal node labelled `inf` or `NaN` is read as support, loses its label, and gets a generated `NODE_` name

### Example data

The input trees of `data/ebola/362` and `data/flu/h3n2/{20,200,500}` carry support values between 0 and 1 with three decimals, the form of FastTree SH-like local support. The datasets do not record how the trees were built. The other example datasets carry no support values.

## Current state in v1

- **Parsing**: `fn graph_from_newick()` ([packages/treetime-io/src/nwk.rs](../../packages/treetime-io/src/nwk.rs)) discards the support value of each node. An internal node whose label was a support value has no name, so `fn assign_node_names()` gives it a `NODE_` name, as v0 does
- **Outputs**: no output writes support: no Auspice `confidence` node attribute or coloring, no `confidence` field in augur node data, no labels or comments in Newick and Nexus
- **Computation**: no step of the core reads support values. They do not affect any inferred value

## Situation in other packages

### TreeTime v0

- **Reading**: Biopython stores the value as `clade.confidence`. No v0 computation reads it
- **Tree files**: the export of `ancestral`, `timetree` and `arg` sets `n.confidence = None` before it writes Newick and Nexus ([packages/legacy/treetime/treetime/CLI_io.py#L160](../../packages/legacy/treetime/treetime/CLI_io.py#L160)), and `mugration` does the same ([packages/legacy/treetime/treetime/wrappers.py#L897](../../packages/legacy/treetime/treetime/wrappers.py#L897)). Only `clock` writes the values, fused into the node names and on the wrong splits after a reroot ([kb/v0-errata/clock-newick-support-fused-into-node-names.md](../v0-errata/clock-newick-support-fused-into-node-names.md))
- **Auspice JSON**: v0 writes a mutation-based substitute under `node_attrs.confidence`, with the coloring `{'title': 'Branch Support', 'type': 'continuous', 'key': 'confidence'}` ([packages/legacy/treetime/treetime/CLI_io.py#L300](../../packages/legacy/treetime/treetime/CLI_io.py#L300), [packages/legacy/treetime/treetime/CLI_io.py#L318-L327](../../packages/legacy/treetime/treetime/CLI_io.py#L318-L327)). For an internal node whose branch carries $m$ substitutions to `A`, `C`, `G` or `T`,

  $$c = 1 - e^{-m}$$

  rounded to three decimals; every leaf gets $c = 1$. v0 writes it only when run with sequence data. The code comment calls it "the bootstrap confidence for iid mutations": when the alignment columns are resampled, the number of the $m$ supporting sites drawn is about Poisson($m$), so the probability that a replicate keeps at least one of them, and with it the split, is $1 - e^{-m}$

### augur

- **`augur refine`** reads the tree with Biopython and writes `clade.confidence` of every node into its node data as `confidence` ([[src](https://github.com/nextstrain/augur/blob/0d287496eed3816f94d674f5a08201273f479487/augur/refine.py#L234)] [[src](https://github.com/nextstrain/augur/blob/0d287496eed3816f94d674f5a08201273f479487/augur/refine.py#L95-L103)])
- **Reroot**: `augur refine --timetree` reroots by default (`--root best`, [[src](https://github.com/nextstrain/augur/blob/0d287496eed3816f94d674f5a08201273f479487/augur/refine.py#L187)]) through TreeTime v0, and outgroup and mid-point rooting also go through Biopython. Biopython's `Tree.root_with_outgroup()` moves branch lengths but not `confidence` ([[src](https://github.com/biopython/biopython/blob/d7e4b8b19399668b09442a5b35765d9186b5f665/Bio/Phylo/BaseTree.py#L873-L878)]), so augur writes the values of rerooted trees on the wrong splits. On the example below, Biopython 1.88 puts 80 above `{A, B, C}`, nothing above `{A, B}`, and 70 above the leaf branch that separates `{A, B, C, D}` from `E`
- **`augur export v2`** treats `confidence` as an ordinary node-data trait, adds an automatic `confidence` coloring and writes `node_attrs.confidence.value` ([[src](https://github.com/nextstrain/augur/blob/0d287496eed3816f94d674f5a08201273f479487/augur/export_v2.py#L361-L362)] [[src](https://github.com/nextstrain/augur/blob/0d287496eed3816f94d674f5a08201273f479487/augur/export_v2.py#L865-L899)])
- **`augur tree`** runs FastTree with `-nosupport`, IQ-TREE without `-B` or `-alrt`, and RAxML with `-f d` by default ([[src](https://github.com/nextstrain/augur/blob/0d287496eed3816f94d674f5a08201273f479487/augur/tree.py#L29-L46)]), so standard Nextstrain builds carry no support values. They appear only in trees built elsewhere or with custom builder arguments
- **Similar names**: `num_date_confidence` (a date interval) and `<trait>_confidence` (state probabilities) are different values. augur writes them as `confidence` sub-fields of their own attributes

### Auspice

Auspice offers a node attribute in its Color By menu only when `meta.colorings` lists it ([[src](https://github.com/nextstrain/auspice/blob/37bf9ce1e3b9a8cbdf1ebfd77bd8d7bdb6dfdc2d/src/components/controls/color-by.js#L210-L211)]). A support value written without a coloring entry cannot be selected in the viewer.

### Phylogenetics tools

- **gotree** stores support on the edge object ([[src](https://github.com/evolbioinfo/gotree/blob/69f84c75d02b5463e872c8513efadd7326fb56b9/tree/edge.go#L20)]), so `Reroot()` only flips edge direction and each value stays with its split. When `UnRoot()` removes a root with two children, the merged edge gets the larger of the two values ([[src](https://github.com/evolbioinfo/gotree/blob/69f84c75d02b5463e872c8513efadd7326fb56b9/tree/tree.go#L1486-L1487)])
- **IQ-TREE** (`-sup`, [[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/utils/tools.cpp#L1630)]) and **RAxML-NG** (`--support`, [[src](https://github.com/amkozlov/raxml-ng/blob/d396351ee1263b704adb135edd6a4a84522dbd5b/src/CommandLineParser.cpp#L1720)]) put support values onto a given tree by matching splits

## Problems to solve before writing support again

- **Placement**: a support value belongs to a split, but the runs edit the topology. `clock` and `timetree` reroot by default and `optimize` reroots with `--reroot`; a reroot inverts the edges between the old and the new root, splits the branch that holds the new root, and removes the old root when it is left with one child (`fn apply_reroot_topology()`, `fn split_edge()`, `fn remove_node_if_trivial()` in [packages/treetime-graph/src/reroot.rs](../../packages/treetime-graph/src/reroot.rs)). `optimize` and `prune` collapse branches, `prune` removes leaves ([packages/treetime/src/prune/prune.rs](../../packages/treetime/src/prune/prune.rs)), and `timetree` resolves polytomies ([packages/treetime/src/timetree/optimization/polytomy/apply.rs](../../packages/treetime/src/timetree/optimization/polytomy/apply.rs)). Values keyed by node, as the Newick reader produces them, end up on the wrong splits:
  - Example: the input root has children `A`, `B` and `X = (C, (D, E)70)80`. After a reroot on the branch to `E`, the node above `{A, B, C}` would show 80, although 80 was measured for `{A, B}` against `{C, D, E}`; its correct value is 70, and the node above `{A, B}` would show nothing instead of 80
  - Measured: with values keyed by node, `treetime optimize --reroot=min-dev` on `data/flu/h3n2/20` places 7 of 14 values on the wrong split, each moved one step along the reroot path
- **Meaning of the Auspice `confidence` attribute**: v0 writes the mutation-based value, augur writes the input support. v1 must choose what the key holds and whether both values are written
- **Visibility in Auspice**: every written support value needs a `meta.colorings` entry
- **Chained runs**: v1 Newick and Nexus outputs carry no support, so a run that reads the tree of an earlier run never sees the input values. `data/ebola/362/pipeline.yaml` passes the Newick output of `optimize` to `ancestral` and `timetree`
- **Label traps**: combined labels such as `80.5/95` become node names, and `inf` or `nan` labels are read as support

## Decided placement rules

Support is stored per edge: the Newick parse keys each value by the edge above its node and drops a value on the root, which has no branch. Each topology edit then follows one rule (approved 2026-10-05):

- An inverted edge keeps its key, so its value stays
- A new root placed inside a branch (`fn split_edge()`): both halves get the branch's value
- Removing a node with one parent and one child (`fn remove_node_if_trivial()`: the old root after a reroot, polytomy cleanup): the merged branch gets the larger of the two values, or none when it ends at a leaf, as gotree does when it removes a root with two children
- Removing a root with one child (`fn remove_stem_root()`): the removed branch's value is dropped
- Collapsing a branch (`fn collapse_edge()` in [packages/treetime/src/optimize/topology/collapse.rs](../../packages/treetime/src/optimize/topology/collapse.rs)): its value is dropped; the child branches keep their keys and values
- A branch created for a new split (polytomy resolution): no value
- Reordering children (`topology_order.apply()`) keeps edge keys, so values stay

## Possible solutions

### How support follows the topology edits of a run

> [!IMPORTANT]
> **Decision required.** The rules above say where a value ends up; the mechanism that applies them is open. The pipelines that edit topology (`clock`, `timetree`, `optimize`, `prune`) do not see support, and the core never computes with it.
>
> - **Edit tracking in the pipelines**: each pipeline takes a per-edge support map as input, updates it wherever it updates branch lengths for a topology edit (beside `fn record_split()` and `fn record_merge()` in [packages/treetime-graph/src/reroot.rs](../../packages/treetime-graph/src/reroot.rs)), and returns it. Small and follows the branch-length pattern; the core carries a value it never computes with, and every edit site must update it
> - **Edit log replayed by the runners**: the pipelines return a log of their topology edits, and the runners replay it on the support map. The core never sees support; every edit, including the polytomy merges inside the iterations of `timetree`, must be recorded exactly
> - **One per-edge record**: one record holds branch length and support, and every edit updates it once. The cleanest data model; it changes every branch-length user in the core
> - **Matching by split after the run**: before the run, record the set of leaves below each input branch with its value; after the run, give each branch of the final tree the value of the same split. This is how IQ-TREE `-sup` and RAxML-NG `--support` place values. No pipeline changes, and the rules above follow from split identity: an inverted branch and both halves of a split branch keep the same split, a collapsed branch's split disappears, a new polytomy branch has a split the input lacks, and two merged input branches map to one split, which takes the larger value. One difference: when a run collapses a split and later builds the same grouping again, matching gives back the input value, while edit tracking gives none. It requires leaf keys that stay stable through every pipeline (reroot keeps node keys, and no pipeline builds a new graph) and costs one leaf bitset per branch. Keying splits by leaf node keys instead of names also avoids the duplicate-name problem

### What Auspice `confidence` holds

> [!IMPORTANT]
> **Decision required.** v0 and augur write different quantities under the same key, and v1 writes neither.
>
> - **Input support**, with a coloring: matches augur (a FastTree value `0.959` appears as `confidence`); a divergence from v0 that needs a `kb/decisions/` entry
> - **The mutation-based value of v0**, with the `Branch Support` coloring: exact v0 parity; trees without input support also get a value (a branch with two substitutions shows $1 - e^{-2} \approx 0.865$)
> - **Both, under separate keys**, each with its own coloring: `confidence` holds input support as in augur, and the mutation-based value gets a new key, which needs a name. Most Nextstrain trees carry no input support, so the mutation-based value is their only support signal; one key for both would change what the number means from input to input

### Support in Newick and Nexus outputs

> [!IMPORTANT]
> **Decision required.** v0 writes no support in its tree files (except the defective `clock` output), and v1 names every internal node. Options: keep writing none; write support as a comment in the annotated styles (`[&support=0.95]` with `--output-nwk-style beast`); or write the value as the label of internal nodes that had no name in the input, which conflicts with the `NODE_` names that the other outputs use to join on internal nodes.

### Label parsing

> [!IMPORTANT]
> **Decision required.** Whether to split combined IQ-TREE labels into several values or keep them as names, and whether `inf` and `nan` labels should be read as names instead of support.

## Validation

- Unit tests for each placement rule on the example above: 80 above `{A, B}` and 70 above `{A, B, C}` after the reroot; a reroot inside a branch copies its value to both halves; removing a root with two children keeps the larger value; a collapse drops the collapsed branch's value; polytomy branches get none
- An oracle independent of the edit path: compare each written value with the input value of the same split, on `data/flu/h3n2/20`, `data/flu/h3n2/200` and `data/ebola/362`, for every command that reroots, collapses or prunes
- Auspice output passes the augur v2 schema with the chosen coloring
- `dev/smoke` cases on the datasets with support values

## References

1. <a id="ref-1"></a>Felsenstein J. 1985. "Confidence Limits on Phylogenies: An Approach Using the Bootstrap." _Evolution_ 39(4):783-791. https://doi.org/10.2307/2408678
2. <a id="ref-2"></a>Minh BQ, Nguyen MAT, von Haeseler A. 2013. "Ultrafast Approximation for Phylogenetic Bootstrap." _Molecular Biology and Evolution_ 30(5):1188-1195. https://doi.org/10.1093/molbev/mst024
3. <a id="ref-3"></a>Anisimova M, Gascuel O. 2006. "Approximate Likelihood-Ratio Test for Branches: A Fast, Accurate, and Powerful Alternative." _Systematic Biology_ 55(4):539-552. https://doi.org/10.1080/10635150600755453
4. <a id="ref-4"></a>Price MN, Dehal PS, Arkin AP. 2010. "FastTree 2 - Approximately Maximum-Likelihood Trees for Large Alignments." _PLoS ONE_ 5(3):e9490. https://doi.org/10.1371/journal.pone.0009490
5. <a id="ref-5"></a>Lemoine F, Domelevo Entfellner JB, Wilkinson E, et al. 2018. "Renewing Felsenstein's phylogenetic bootstrap in the era of big data." _Nature_ 556:452-456. https://doi.org/10.1038/s41586-018-0043-0
6. <a id="ref-6"></a>Czech L, Huerta-Cepas J, Stamatakis A. 2017. "A Critical Review on the Use of Support Values in Tree Viewers and Bioinformatics Toolkits." _Molecular Biology and Evolution_ 34(6):1535-1542. https://doi.org/10.1093/molbev/msx055
