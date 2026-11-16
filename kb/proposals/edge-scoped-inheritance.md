# Edge-scoped inheritance: per-site local trees on one graph

Segmented viruses reassort: a progeny genome can take one segment from parent A and another from parent B. Recombining genomes go further, so that a single segment takes the sites left of a crossover point from one parent and the sites right of it from the other. In both cases every site still has exactly one ancestor at each node, so every site has its own <a id="gloss-use-1"></a>local tree <sup>[1](#gloss-1)</sup>, but different sites have different local trees.

This proposal represents such histories on the existing TreeTime graph. Each edge records which sites are inherited along it, and the partition passes run over the whole graph, sending each site only along the edges that carry it. Trees, reassortment, and crossover recombination then differ only in the edge scopes. Every local tree is used through the scopes; no subtree or forest of the graph is built.

## Current state

- **Graph**: the graph and its schedule already accept directed acyclic graphs (DAGs) with several parents per node and several roots. The schedule rejects cycles and orders every node after all its parents (preorder) or children (postorder), and `scan_down` starts one task per root ([kb/algo/graph.md](../algo/graph.md))
- **Partition passes**: Fitch and the dense, discrete, and sparse marginal passes call `require_tree()` (one root, at most one parent per node) at their start (`fn require_tree()` [packages/treetime-graph/src/schedule.rs#L54](../../packages/treetime-graph/src/schedule.rs#L54)). Their forward steps take the message of the first parent only (`parents.first()`, for example [packages/treetime/src/partition/marginal/shared/pass.rs#L281](../../packages/treetime/src/partition/marginal/shared/pass.rs#L281) and [packages/treetime/src/partition/fitch/passes.rs#L209](../../packages/treetime/src/partition/fitch/passes.rs#L209)), and their backward steps treat a node without parents as the root of all sites ([packages/treetime/src/partition/marginal/shared/pass.rs#L193](../../packages/treetime/src/partition/marginal/shared/pass.rs#L193))
- **Earlier transmission ranges**: the sparse edge state had an optional `transmission: Option<Vec<(usize, usize)>>` of site ranges, meant to record what an edge inherits. No code ever assigned it, and its two readers disagreed: the Fitch backward step skipped a child at a position inside the ranges, while `propagate_raw` skipped a position outside them. Commit `cab785f66` (`refactor(sparse): remove the never-assigned edge transmission ranges`) removed it
- **Input**: one tree per run. Segments with different trees need one run each ([kb/issues/N-io-multi-segment-genome-input.md](../issues/N-io-multi-segment-genome-input.md)); the `arg` command of v0 is not ported ([kb/features/arg.md](../features/arg.md))

### v0 masks

v0 `treetime arg` analyses one focal segment tree at a time with a concatenated alignment of all segments. It assigns each node to a <a id="gloss-use-2"></a>maximum compatible clade (MCC) <sup>[2](#gloss-2)</sup> of TreeKnit and gives every branch a mask: the whole alignment when parent and child lie in the same MCC, the focal segment only otherwise (`fn setup_arg()` `packages/legacy/treetime/treetime/arg.py:128-133`). The marginal backward pass multiplies the log-messages of masked sites by zero, so a masked site sends a uniform message across that branch (`packages/legacy/treetime/treetime/treeanc.py:862-871`), and the forward pass replaces the parent message at masked sites by one (`packages/legacy/treetime/treetime/treeanc.py:909-916`).

In the terms of this proposal, v0 uses the focal tree as the graph, the focal segment has full scope on every edge, and each other segment has scope only on the edges inside an MCC. The local trees of the other segments are cut at every MCC root, and the cut subtrees have no root prior. The masks also enter branch-length optimization and the timetree through the pattern multiplicity (`packages/legacy/treetime/treetime/clock_tree.py:363`).

## Proposed model

### Edge scope

An edge scope says which sites of which partition an edge carries. It lives in a value map keyed by edge (`EdgeMap<EdgeScope>`), next to the branch lengths, so the graph stays topology-only ([kb/decisions/graph-topology-only-value-pipeline.md](../decisions/graph-topology-only-value-pipeline.md)).

- **Per partition**: a set of half-open site ranges in the coordinates of the partition. Empty means the edge does not carry the partition, full means it carries every site
- **Tree input**: every edge has full scope in every partition. This is the default when no scope map is given, so tree inputs run exactly as today
- **Reassortment**: every range set is empty or full: an edge carries a segment entirely or not at all
- **Crossover recombination**: the two parent edges of a recombinant node split a partition, for example `[0, x)` on one and `[x, L)` on the other

### Invariants

For every node `v`, partition `p`, and site `i` of `p`:

- **One parent per site**: at most one parent edge of `v` carries `i`, so the scopes of the parent edges of `v` are disjoint
- **No inheritance without descendants**: a parent edge of `v` carries `i` only if `v` is a leaf with data for `p` or some child edge of `v` carries `i`

The graph must be acyclic, which the schedule checks. These two invariants replace the tree check of the partition passes. Each node then holds a set of sites per partition: those carried by its child edges, or all sites of `p` at a leaf with data for `p`. A site that a node holds but receives from no parent is a <a id="gloss-use-3"></a>root site <sup>[3](#gloss-3)</sup> of that node. A node may be a root for some sites and an inner node for others, and a partition may have any number of root nodes. Roots follow from the scopes and need no search.

### Passes

The passes keep the graph schedule and its two scans; only the node steps read the scopes.

- **Backward (`scan_up`)**: a node combines the message of each child edge over that edge's scope and treats sites outside it as giving no information (factor one). It sends each parent edge its message restricted to that edge's scope; the scans already deliver one message per parent slot. At its root sites it multiplies by the equilibrium frequencies and adds the normalization to the log-likelihood, as the current root step does for all sites
- **Forward (`scan_down`)**: a node builds its message from above site by site, taking each site from the parent edge that carries it, and uses the equilibrium frequencies at its root sites
- **Fitch**: the same rules with state sets in place of profiles. A site outside a child's scope is skipped like a site where the child has no data (`range_contains(&child.non_char, pos)` [packages/treetime/src/partition/fitch/sub.rs#L36](../../packages/treetime/src/partition/fitch/sub.rs#L36))
- **Sparse**: a site outside a child edge's scope is a site without data below that child, which the sparse backward pass already handles through the no-data ranges of each child (`child_no_data` [packages/treetime/src/partition/storage/sparse.rs#L142](../../packages/treetime/src/partition/storage/sparse.rs#L142)). The forward pass needs the per-site choice of parent. Whole-partition scopes need no change to the variable-site maps, because an edge then carries the whole partition or is skipped
- **Discrete traits (mugration)**: a trait is one site. A trait partition takes the scope of one segment, or full scope, and is not split
- **Log-likelihood**: the sum over all root sites of all nodes. On a tree this is the root term of today

For every site the edges that carry it form a tree, so marginal and Fitch reconstruction stay exact. For reassortment with whole-segment scopes, a partition gives the same states, mutations, and log-likelihood as a run on its segment tree alone; this is the main validation.

### Branch lengths

An edge has one branch length, shared by all partitions it carries, as on a tree. Partitions that evolve at different speeds use the per-partition rate (`enum InferredRate` [packages/treetime/src/partition/inferred_rate.rs](../../packages/treetime/src/partition/inferred_rate.rs)). The branch-length objective of an edge sums the contributions of the partitions and sites the edge carries ([packages/treetime/src/partition/optimize/contribution.rs](../../packages/treetime/src/partition/optimize/contribution.rs)), so optimization needs the scope but no new structure. One length per edge is also what the timetree needs, where an edge is a time interval.

### Time trees

A node has one time, shared by all its sites. A reassortant node with two older parents closes an undirected cycle, so the time inference on the graph is no longer a tree problem: two passes of message passing give exact marginals on a tree, but on a graph with undirected cycles they become loopy belief propagation, which is approximate <a id="cite-1"></a>[Koller and Friedman 2009](https://mitpress.mit.edu/9780262013192/probabilistic-graphical-models/) [[1](#ref-1)]. Sequence inference is not affected, because each site sees a tree.

> [!IMPORTANT]
> **Decision required.** Time inference on a graph with reassortment needs one of:
>
> - Loopy message passing on the time graph, iterated to a fixed point. Reassortment nodes are few compared to all nodes, so the cycles are few, but convergence and calibration need checking
> - Conditioning on one parent edge per reassortant node in turn, as v0 does with its focal trees, and combining the results
> - Exact treatment of the cycles, for example by grouping the nodes of each cycle (junction tree), with a cost that grows with the size of the cycles
>
> The choice also settles which timetree parameters the segments share (clock rate, coalescent model), the open question of [kb/issues/N-io-multi-segment-genome-input.md](../issues/N-io-multi-segment-genome-input.md).

Decided: loopy message passing, node-local. The backward and forward passes run unchanged at nodes with one parent; a node with several parents sends each parent its subtree (optionally combined, from the second sweep on, with the messages of its other parents), and multiplies the messages of all parents into its posterior. Sweeps repeat until the node times settle ([kb/algo/timetree.md](../algo/timetree.md#graphs-with-several-parents)). The segments share one clock rate, given with `--clock-rate`.

### Topology operations

Rerooting, polytomy resolution, and edge contraction change edges, and the scopes must follow. Contracting an edge `X -> Y` gives each child edge `Y -> c` the new source `X`, which is valid for the sites of `Y -> c` that `Y` inherits from `X`, but not for sites that `Y` inherits from another parent. Rerooting inverts edges, which has no meaning at a reassortant node. A first version runs these operations only when all scopes are full (a tree) and reports an error otherwise.

## Relation to other representations

- **Subtree per partition**: building one tree per partition from the scopes, with node and edge keys kept, would let the current tree passes run unchanged. It does not extend to crossover recombination, where a partition has a different tree in each interval, and it splits the shared node times of the timetree across partitions. This proposal keeps one graph and one pass per partition instead
- **v0 masks**: the special case above; the edge scope replaces the per-branch mask, and the graph replaces the focal tree
- **Tree sequences**: tskit stores an ARG as edges with a left and right genome coordinate, a parent, and a child <a id="cite-2"></a>[Kelleher et al. 2018](https://doi.org/10.1371/journal.pcbi.1006581) [[2](#ref-2)]. An edge scope with several ranges is the same information with the ranges grouped by graph edge; the general ARG representation of <a id="cite-3"></a>[Wong et al. 2024](https://doi.org/10.1093/genetics/iyae100) [[3](#ref-3)] likewise attaches genome intervals to edges

## Input

- **Named trees**: one Newick tree per segment, with named internal nodes. A name that occurs in several trees is one graph node, and a parent-child pair that occurs in several trees is one graph edge whose scope holds all their segments. A node that a segment only passes through belongs in that segment's tree as a node with one child, so that an edge shared by several segments stays one edge; the merge inserts such nodes where topology determines them. A shared edge gets the mean of the branch lengths its trees give, since separately estimated segment trees need not agree
- **Annotated network Newick**: an extended Newick network whose edges carry segment sets, such as `[&segments={0,1}]` of CoalRe, holds the same information in one file; the tests of the Newick dialects already read such a string (`COALRE` in [packages/app-commands/src/commands/shared/__tests__/test_tree_input.rs](../../packages/app-commands/src/commands/shared/__tests__/test_tree_input.rs))
- **TreeKnit**: TreeKnit writes one tree per segment with names that are local to each tree, and the MCCs of every pair of segments (`data/flu/treeknit/MCCs.json`, four segments and six pairs). Turning these into named trees needs a rule for MCCs of different pairs that disagree, and the nodes that segments pass through must be inserted. A prototype matches nodes as MRCAs of the same MCC leaves and runs on the H3N2 data; its results and the information that TreeKnit does not provide are in [kb/reports/2026-11-15-treeknit-output-to-named-segment-trees.md](../reports/2026-11-15-treeknit-output-to-named-segment-trees.md)
- **Toy data**: [data/toy/reassortment/](../../data/toy/reassortment/) holds a hand-made graph with three segments, two reassortant nodes (one with a different parent for each segment), two graph roots, and nodes that some segments pass through, with simulated alignments and the true ancestral sequences

## Open design axes

- **Scope storage** (implemented in `treetime-graph`, [kb/algo/graph.md](../algo/graph.md#edge-scopes)): one `SiteScope` per segment and edge, `Empty`, `Full`, or a few ranges stored inline, so reassortment and a few crossovers per segment do not allocate. Many breakpoints per segment, as in human ARGs, are out of scope
- **Coordinates** (implemented): ranges are in the coordinates of a segment, which is one segment tree with its alignment. An alignment split into several partitions, such as genes, maps the segment ranges to the columns of each partition
- **Command interface** (first version implemented): `ancestral --segment-trees ... --segment-alignments ...`, paired in order, separate from `--tree` and `--alignment`; one interface with them and the YAML configuration ([kb/proposals/config-file-multi-partition.md](config-file-multi-partition.md)) later
- **Outputs**: the graph with scopes as graph JSON; per segment, the tree that the segment's scopes trace, in the tree formats

## Work items

- Read named segment trees into one graph with edge scopes: [kb/issues/N-io-multi-segment-genome-input.md](../issues/N-io-multi-segment-genome-input.md)
- Scope-aware Fitch and marginal passes: [kb/issues/M-partition-passes-crossover-scopes-and-discrete-traits.md](../issues/M-partition-passes-crossover-scopes-and-discrete-traits.md)
- Branch-length optimization on graphs of segment trees: [kb/issues/M-optimize-segment-graphs-topology-fixed.md](../issues/M-optimize-segment-graphs-topology-fixed.md)
- Time inference on graphs with reassortment: [kb/issues/M-timetree-segment-graphs-fixed-clock-only.md](../issues/M-timetree-segment-graphs-fixed-clock-only.md)

## Glossary

1. <a id="gloss-1"></a> **Local tree.** The tree that the ancestry of one site traces through a graph with reassortment or recombination. [↩](#gloss-use-1)
2. <a id="gloss-2"></a> **Maximum compatible clade (MCC).** In TreeKnit, a maximal set of leaves whose subtree topology two segment trees share; reassortment events lie on the branches above MCC roots. [↩](#gloss-use-2)
3. <a id="gloss-3"></a> **Root site.** A site that a node holds (through its children or its own data) but inherits from none of its parents. [↩](#gloss-use-3)

## References

1. <a id="ref-1"></a> Koller, Daphne, and Nir Friedman. 2009. _Probabilistic Graphical Models: Principles and Techniques._ MIT Press. ISBN 978-0-262-01319-2. https://mitpress.mit.edu/9780262013192/probabilistic-graphical-models/ [↩](#cite-1)
2. <a id="ref-2"></a> Kelleher, Jerome, Kevin R. Thornton, Jaime Ashander, and Peter L. Ralph. 2018. "Efficient pedigree recording for fast population genetics simulation." _PLOS Computational Biology_ 14:e1006581. https://doi.org/10.1371/journal.pcbi.1006581 [↩](#cite-2)
3. <a id="ref-3"></a> Wong, Yan, Anastasia Ignatieva, Jere Koskela, Gregor Gorjanc, Anthony W. Wohns, and Jerome Kelleher. 2024. "A general and efficient representation of ancestral recombination graphs." _Genetics_ 228:iyae100. https://doi.org/10.1093/genetics/iyae100 [↩](#cite-3)
