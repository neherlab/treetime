# Graph-based phylogenetic representation

v1 replaces BioPython's tree-only data model with a directed graph capable of representing both phylogenetic trees and phylogenetic networks. `Graph<N, E, D>` (`#Graph`) in [packages/treetime-graph/src/graph.rs#L32-L43](../../packages/treetime-graph/src/graph.rs#L32-L43) replaces `Phylo.BaseTree.Tree` wrapped by `TreeAnc` (`#TreeAnc`) in [packages/legacy/treetime/treetime/treeanc.py#L50-L514](../../packages/legacy/treetime/treetime/treeanc.py#L50-L514).

## Background: phylogenetic trees

A phylogenetic tree is a graph-theoretical representation of evolutionary relationships among taxa (Qi & Schicho, 2020 [[1]](#ref-1)). Nodes represent either observed samples (leaves) or inferred ancestors (internal nodes). Edges represent evolutionary lineages with associated branch lengths measured in expected substitutions per site.

Five mathematically equivalent representations of a phylogenetic tree exist: edge-labeled graphs, split systems, hierarchical cluster systems, tree metrics, and systems of nested parentheses (Qi & Schicho, 2020 [[1]](#ref-1)). This equivalence means any representation can be converted to any other without loss, and the choice is a matter of algorithmic convenience.

Most phylogenetic software uses the Newick format for tree serialization: `(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);`. The Newick Standard was adopted June 26, 1986 by an informal committee meeting at Newick's restaurant in Dover, New Hampshire (Felsenstein [[2]](#ref-2)). The format generalizes an earlier notation developed by Christopher Meacham in 1984 for PHYLIP (Felsenstein [[2]](#ref-2)). It traces mathematically to Arthur Cayley's 1857 observation of the correspondence between trees and nested parentheses (Felsenstein [[2]](#ref-2)).

The Newick format encodes parent-child relationships implicitly through nesting. The standard data structure mirrors this: each node stores a list of children, and branch length is stored on the child node (Felsenstein [[2]](#ref-2)). BioPython's `Bio.Phylo.BaseTree.Clade` follows this convention (Talevich et al., 2012 [[3]](#ref-3)).

This child-linked representation has two limitations for algorithms that traverse toward the root. First, accessing a node's parent requires either a separate lookup or an explicit back-pointer. Second, branch data conceptually belongs to the connection between nodes, not to either node individually, but the format forces it onto the child.

## The tree assumption and its limitations

Standard phylogenetic trees assume strictly divergent (bifurcating) evolution: each lineage has exactly one ancestor, and lineages never merge after splitting. This assumption holds for many macroevolutionary scenarios but breaks down in several biologically important cases.

### Reticulate evolution

Reticulate evolution is the partial merging of ancestor lineages, creating relationships that cannot be depicted as a bifurcating tree (Huson et al., 2011 [[4]](#ref-4)). It occurs across multiple biological scales.

**Hybridization and polyploidy.** Roughly 25% of plant species and 10% of animal species form hybrids with at least one other species (Gontier, 2015 [[5]](#ref-5)). Allopolyploidy (whole-genome duplication following hybridization between different species) is a major speciation mechanism in plants. Bread wheat is hexaploid, arising from three successive hybridization events. Cotton, canola, and peanut are allopolyploids. An allopolyploid species has two distinct parental lineages: its homoeologous chromosomes trace back to different ancestral species, creating gene tree discordance that no single bifurcating tree can represent (Gontier, 2015 [[5]](#ref-5)).

**Admixture and introgression.** Introgression is the transfer of genetic material from one species into another through repeated backcrossing after an initial hybridization event. Unlike F1 hybridization (which produces a single reticulation node), introgression transfers alleles gradually over many generations, creating a diffuse signal across the genome. The Neanderthal genome project (Green et al., 2010 [[31]](#ref-31)) demonstrated that 1-4% of the genome of non-African modern humans derives from Neanderthal introgression. The D-statistic (ABBA-BABA test) detects introgression by comparing site patterns across four taxa: a significant excess of shared derived alleles between non-sister taxa indicates gene flow beyond what ILS alone would produce. Introgression is common in plant evolution, fish radiations (cichlids), and canid hybridization (wolves, coyotes, dogs).

**Horizontal gene transfer.** In prokaryotes, genetic material moves between organisms through three mechanisms: transformation (uptake of free DNA), transduction (phage-mediated transfer), and conjugation (direct cell-to-cell transfer via pilus) (Ge et al., 2005 [[6]](#ref-6)). About 2% of core genes per prokaryotic genome show evidence of horizontal transfer (Ge et al., 2005 [[6]](#ref-6)). This rate is low enough that a species tree remains meaningful but high enough that gene-level phylogenies frequently conflict with it. Horizontal gene transfer is the primary mechanism for spreading antibiotic resistance among bacteria.

**Recombination in viruses.** RNA viruses undergo frequent genetic recombination via copy-choice (template switching), where the RNA-dependent RNA polymerase switches template during replication and produces chimeric progeny (Wikipedia [[7]](#ref-7)). SARS-CoV-2 variants designated with an "X" prefix (XBB, XEC, XFG) are recombinant lineages derived from co-infection with distinct parent lineages (Wikipedia [[7]](#ref-7)). Different genomic regions of a recombinant virus yield different tree topologies when analyzed independently.

**Reassortment in segmented viruses.** Influenza A has eight genome segments. When two strains co-infect a cell, segments from different parents can be packaged together, producing reassortant progeny with mixed-origin genomes (Wikipedia [[8]](#ref-8)). All major influenza pandemics since the early 1900s involved reassortment: the 1957 Asian Flu acquired HA and NA segments from avian influenza, the 1968 Hong Kong Flu reassorted a new HA from an avian strain, and the 2009 pandemic strain combined segments from four distinct lineages (avian, two swine, and human) (Wikipedia [[8]](#ref-8)). Pigs serve as "mixing vessels" because their respiratory epithelium expresses both avian-type and human-type sialic acid receptors, allowing co-infection by viruses from multiple host species.

### Incomplete lineage sorting

Even without reticulation, gene trees can disagree with species trees through incomplete lineage sorting (ILS) (Wikipedia [[9]](#ref-9)). When ancestral populations are polymorphic at speciation events and the inter-speciation interval is short relative to the effective population size, different loci can fix different alleles in descendant species. For a rooted three-taxon species tree with internal branch length T (in coalescent units), the probability that a gene tree matches the species tree is 1 - (2/3)exp(-T) (Wikipedia [[9]](#ref-9)). About 23% of DNA alignments in great apes do not support the known human-chimpanzee sister relationship (Wikipedia [[9]](#ref-9)). ILS does not create reticulate topology in the species tree itself but produces the same data-level pattern of conflicting topologies across loci.

### Gene duplication and loss

Gene duplication creates paralogous copies within a genome (Wikipedia [[32]](#ref-32)). If a duplication occurs before a speciation event and one copy is subsequently lost in some lineages (differential gene loss), orthology inference can mistake a paralog for an ortholog. The resulting gene tree reflects the duplication history rather than the species history, producing a topology that disagrees with the species tree. Unlike HGT and hybridization, gene duplication does not create reticulation in the species tree itself, but it generates the same data-level pattern of conflicting gene trees. Phylogenetic reconciliation methods distinguish duplication/loss events from speciation events by mapping gene trees into species trees and minimizing the total number of inferred events.

## Phylogenetic networks

A phylogenetic network is a graph used to represent evolutionary relationships that include reticulation events (Huson et al., 2011 [[4]](#ref-4)). Trees are a special case of networks where every node has at most one parent. In a rooted phylogenetic network:

- Tree nodes have exactly one parent (divergence events).
- Hybrid nodes (reticulate nodes) have two or more parents (hybridization, HGT, recombination events).
- Leaves are bijectively labeled by taxa.
- The graph is a directed acyclic graph (DAG) (Huson et al., 2011 [[4]](#ref-4)).

### Network classes

Network classes form a containment hierarchy, each constraining the allowed reticulation patterns (Huson & Scornavacca, 2011 [[10]](#ref-10)):

- Trees (strictest): no reticulate nodes. Each node has at most one parent.
- Galled trees: each reticulate node belongs to a single "gall" (biconnected cycle), and galls share no edges.
- Tree-child networks: every non-leaf node has at least one child that is a tree node (not a reticulate node). Every internal node has a "tree path" to a leaf.
- Normal networks: tree-child with no redundant arcs.
- Reticulation-visible networks: for every reticulate node, there exists a tree node below it that is not below any other reticulate node.
- Tree-based networks: obtainable by adding arcs between arcs of some base tree (Pons et al., 2018 [[11]](#ref-11)). Contains tree-child and reticulation-visible as subclasses.
- Level-k networks: the maximum number of reticulate nodes in any biconnected component is at most k (Huson & Scornavacca, 2011 [[10]](#ref-10)). Level-0 networks are trees. Level-1 networks are close to galled trees.

The treewidth of level-k networks is bounded by (k+3)/2 (Markin et al., 2024 [[12]](#ref-12)), making many NP-hard problems polynomial for constant k.

### Displayed trees

A rooted phylogenetic network displays a set of trees. For each reticulate node, selecting exactly one parent edge (choosing which lineage contributed to the character at that site) produces a tree embedded within the network. The set of all such selections yields all displayed trees. This concept is central to softwired parsimony, where the optimal score is the minimum parsimony score over all displayed trees (Fischer et al., 2015 [[13]](#ref-13)).

### Splits and split networks

A split is a bipartition of a set of taxa (Wikipedia [[14]](#ref-14)). Each edge of an unrooted phylogenetic tree corresponds to one split. Two splits are compatible if they can coexist in a single tree. When a dataset contains incompatible splits (conflicting phylogenetic signal from recombination, hybridization, or systematic error), a split network represents the conflict as arrays of parallel edges rather than forcing a single tree topology (Wikipedia [[14]](#ref-14)). Split networks are unrooted and primarily visualize ambiguity in data rather than committing to a specific evolutionary scenario. NeighborNet (Bryant & Moulton, 2004 [[33]](#ref-33)) is the most widely used algorithm for constructing split networks from distance matrices. It produces overlapping (non-hierarchical) clusters represented as a planar splits graph. When the distance matrix is tree-like (all splits compatible), NeighborNet returns a tree identical to neighbor joining. When splits conflict, the graph contains parallelogram-shaped regions that visualize the conflicting signal. NeighborNet has been applied in virology, horticulture, comparative linguistics, and archaeology.

### Ancestral recombination graphs

An ancestral recombination graph (ARG) represents the genealogical history of DNA sequences accounting for both coalescence and recombination (Wikipedia [[15]](#ref-15)). Coalescence nodes merge two child lineages into one parent (standard coalescent). Recombination nodes split one child lineage into two parent lineages, each carrying a different segment of the genome. The breakpoint defines which genomic segment follows which ancestral path. Along a chromosome, the ARG induces a sequence of "local trees" (marginal genealogies), each applying to a contiguous segment between recombination breakpoints. Adjacent local trees differ by one subtree-prune-and-regraft (SPR) operation (Wikipedia [[15]](#ref-15)). Practical ARG infrastructure includes msprime for simulation under demographic models, tskit for efficient storage and manipulation of ARGs at population scale (Kelleher et al., 2018 [[35]](#ref-35)), and Moonshine.jl for model-based ARG inference from thousands of haplotypes (Fournier & Larribe, 2025 [[36]](#ref-36)).

## Algorithms on phylogenetic networks

### Parsimony

Two parsimony frameworks exist for networks (Fischer et al., 2015 [[13]](#ref-13)):

**Hardwired parsimony** treats the network as a directed graph and counts state changes along all edges, including reticulation edges. For binary characters, this is polynomial-time solvable via reduction to Multicut. For more states, it is NP-hard but fixed-parameter tractable in the parsimony score (Fischer et al., 2015 [[13]](#ref-13)).

**Softwired parsimony** selects one parent edge per reticulate node per site, extracting a displayed tree, and minimizes the parsimony score over all displayed trees. This is NP-hard even for binary characters but fixed-parameter tractable in the network level (Fischer et al., 2015 [[13]](#ref-13)). Fitch's parsimony algorithm extends to binary tree-child networks with a guaranteed 2-approximation factor for softwired parsimony, interpreting Fitch's algorithm as a primal-dual algorithm (Frohn & Kelk, 2024 [[16]](#ref-16)). The gap between a displayed tree's parsimony score and the optimal softwired score is bounded by (k+1) times the softwired score for level-k networks, independent of alignment length (Doecker et al., 2024 [[17]](#ref-17)).

### Likelihood and Bayesian inference

Maximum likelihood on networks is less developed than parsimony. The multispecies network coalescent (MSNC) model accounts for both hybridization and incomplete lineage sorting (Elworth et al., 2019 [[18]](#ref-18)). PhyloNet implements ML under this model. SNaQ uses quartet concordance factors under the coalescent model extended with hybridization, providing a scalable pseudolikelihood approach (Solis-Lemus & Ane, 2016 [[19]](#ref-19)). Bayesian methods (SpeciesNetwork/BEAST 2) jointly infer network topology, divergence times, and inheritance probabilities from multilocus sequence data, though convergence difficulties arise on larger datasets (Zhang et al., 2018 [[20]](#ref-20)).

### Molecular clock inference

Extending molecular clock inference to networks is substantially harder than on trees because reticulations create multiple evolutionary paths between ancestral and extant taxa. Edge length calibration from genetic distances is possible for known network topologies via least-squares decomposition into invariant and variable components (Xu & Ane, 2024 [[21]](#ref-21)). Convergence-divergence models retain a single principal tree but allow gene flow over arbitrary time frames, providing an alternative to explicit network models when reticulation is best modeled as continuous rather than instantaneous (Mitchell & Holland, 2025 [[22]](#ref-22)). Fast approximate ML clock inference (the approach used by TreeTime on trees) has not been extended to networks.

### Identifiability

Tree-child networks are identifiable under a probabilistic recombination-mutation model (Francis & Moulton, 2018 [[23]](#ref-23)). Galled tree-child networks (beyond level-1) are identifiable from quartet concordance factor data (Allman et al., 2025 [[24]](#ref-24)), the strongest identifiability result for networks to date. Binary tree-child networks and recoverable binary level-2 networks are uniquely determined by their trinets (induced 3-leaf subnetworks) (van Iersel & Moulton, 2014 [[25]](#ref-25)). Level-k networks can be reconstructed from dense triplet sets in O(|T|^{k+1}) time (To & Habib, 2009 [[26]](#ref-26)).

## Serialization: Extended Newick

Extended Newick generalizes standard Newick to encode phylogenetic networks (Cardona et al., 2008 [[27]](#ref-27)). Reticulate nodes are duplicated in the string and annotated with `#` followed by an optional type string and a mandatory integer identifier:

```
(A,B,((C,(Y)x#H1)c,(x#H1,D)d)e)f;
```

Here `x#H1` appears twice. A parser joins these into a single node with two incoming edges, producing a DAG rather than a tree. The type string encodes the reticulation event: `H` for hybridization, `LGT` for lateral gene transfer, `R` for recombination (Cardona et al., 2008 [[27]](#ref-27)). Extended Newick is backward-compatible with standard Newick: a legacy parser treats `x#H1` as an oddly named leaf, producing a valid tree.

The `##` (double hash) convention marks the "acceptor" edge among a reticulate node's parents, distinguishing it from "transfer" edges (Cardona et al., 2008 [[27]](#ref-27)). NHX (New Hampshire X) and Rich Newick are additional extensions adding key-value annotations and rooted/unrooted markers.

v1 does not parse Extended Newick. The Newick reader `nwk_read()` in [packages/treetime-io/src/nwk.rs#L45-L103](../../packages/treetime-io/src/nwk.rs#L45-L103) delegates to the `bio::io::newick` crate, which supports standard Newick only. The Newick writer rejects multiple roots with `unimplemented!()`. Extended Newick parsing and multi-root serialization are prerequisites for loading phylogenetic networks from files.

## v0: BioPython Clade with monkey-patched attributes

v0 represents trees using BioPython's `Bio.Phylo` module (Talevich et al., 2012 [[3]](#ref-3)). The tree is a `Phylo.BaseTree.Tree` object. Each node is a `Bio.Phylo.BaseTree.Clade` instance storing children in a `.clades` list attribute (BioPython docs [[28]](#ref-28)).

BioPython's data model provides no explicit parent attribute (BioPython docs [[28]](#ref-28)). To find a parent, one must traverse from the root using `get_path()`. v0 adds parent pointers by monkey-patching: `_prepare_nodes()` (`#_prepare_nodes`) in [packages/legacy/treetime/treetime/treeanc.py#L461-L493](../../packages/legacy/treetime/treetime/treeanc.py#L461-L493) walks the tree in preorder and sets `c.up = clade` on each child. Every node also receives a back-reference to the owning `TreeAnc` instance via `c.tt = self`.

Algorithm-specific data is attached directly to Clade objects as dynamic attributes. Three properties are patched onto the `Clade` class at module level ([packages/legacy/treetime/treetime/treeanc.py#L45-L47](../../packages/legacy/treetime/treetime/treeanc.py#L45-L47)):

- `Clade.sequence` - full sequence array
- `Clade.cseq` - compressed sequence
- `Clade.mutations` - mutations relative to parent

Other attributes are set during algorithm passes without declaration: `node.branch_length_interpolator`, `node.date_constraint`, `node.time_before_present`, `node.marginal_subtree_LH`, `node.joint_Lx`, `node.bad_branch`. Missing attributes cause runtime `AttributeError`.

Traversal uses BioPython methods (BioPython docs [[28]](#ref-28)):

- `tree.find_clades(order='postorder')` for leaves-to-root (depth-first, children before parent)
- `tree.find_clades(order='preorder')` for root-to-leaves (depth-first, parent before children)
- `tree.get_terminals()` for leaf nodes
- `tree.get_nonterminals()` for internal nodes

All traversals are single-threaded and sequential.

Bio.Phylo assumes strictly tree-shaped data (BioPython docs [[28]](#ref-28)). There is no mechanism for representing nodes with multiple parents, no edge entity separate from nodes, and no support for Extended Newick or phylogenetic network serialization.

## v1: generic directed graph with typed payloads

v1 uses `Graph<N, E, D>` defined in the `treetime-graph` crate. The struct is parameterized over three types:

- `N: GraphNode` - node payload containing algorithm-specific data
- `E: GraphEdge` - edge payload containing branch data
- `D: Sync + Send` - graph-level shared data

### Storage

Nodes are stored in `Vec<Option<Arc<RwLock<Node<N>>>>>` indexed by `GraphNodeKey(usize)`. Removed nodes leave `None` tombstones; keys are never reassigned. Each `Node<N>` (`#Node`) in [packages/treetime-graph/src/node.rs#L73-L79](../../packages/treetime-graph/src/node.rs#L73-L79) tracks both inbound and outbound edge keys:

```rust
pub struct Node<N: GraphNode> {
  key: GraphNodeKey,
  data: Arc<RwLock<N>>,
  outbound_edges: Vec<GraphEdgeKey>,
  inbound_edges: Vec<GraphEdgeKey>,
  is_visited: AtomicBool,
}
```

Parent access uses `inbound_edges`: `graph.parents_of(node)` (`#parents_of`) in [packages/treetime-graph/src/graph.rs#L86-L97](../../packages/treetime-graph/src/graph.rs#L86-L97) retrieves parent nodes by following inbound edge keys to their source nodes. Child access uses `outbound_edges` similarly via `graph.children_of(node)` (`#children_of`) in [packages/treetime-graph/src/graph.rs#L143-L154](../../packages/treetime-graph/src/graph.rs#L143-L154).

Edges are stored in `Vec<Option<Arc<RwLock<Edge<E>>>>>`. Each `Edge<E>` (`#Edge`) in [packages/treetime-graph/src/edge.rs#L71-L83](../../packages/treetime-graph/src/edge.rs#L71-L83) stores source and target node keys plus a payload:

```rust
pub struct Edge<E: GraphEdge> {
  key: GraphEdgeKey,
  source: GraphNodeKey,
  target: GraphNodeKey,
  data: Arc<RwLock<E>>,
}
```

The graph stores root and leaf key vectors (`roots: Vec<GraphNodeKey>`, `leaves: Vec<GraphNodeKey>`) recomputed by `build()` (`#build`) in [packages/treetime-graph/src/graph_ops.rs#L122-L131](../../packages/treetime-graph/src/graph_ops.rs#L122-L131). A root is any node with no inbound edges. A leaf is any node with no outbound edges.

### Typed payloads

Node and edge payloads are concrete types per algorithm stage. `NodeAncestral` (`#NodeAncestral`) and `EdgeAncestral` (`#EdgeAncestral`) in [packages/treetime/src/partition/payload/ancestral.rs#L14-L72](../../packages/treetime/src/partition/payload/ancestral.rs#L14-L72) provide minimal fields for ancestral reconstruction. `NodeTimetree` (`#NodeTimetree`) and `EdgeTimetree` (`#EdgeTimetree`) in [packages/treetime/src/partition/payload/timetree.rs#L18-L158](../../packages/treetime/src/partition/payload/timetree.rs#L18-L158) extend these with time-related fields.

Traits define typed access to domain-specific fields:

- `Named`, `Described` - node identification ([packages/treetime-graph/src/node.rs#L11-L20](../../packages/treetime-graph/src/node.rs#L11-L20))
- `Divergence` - evolutionary distance from root ([packages/treetime-graph/src/node.rs#L32-L35](../../packages/treetime-graph/src/node.rs#L32-L35))
- `TimeConstraint<T>` - date constraints and bad branch flag ([packages/treetime-graph/src/node.rs#L44-L49](../../packages/treetime-graph/src/node.rs#L44-L49))
- `HasBranchLength` - branch length on edges ([packages/treetime-graph/src/edge.rs#L13-L16](../../packages/treetime-graph/src/edge.rs#L13-L16))
- `ClockMessages<T>` - message-passing fields for clock inference ([packages/treetime-graph/src/edge.rs#L20-L27](../../packages/treetime-graph/src/edge.rs#L20-L27))

Functions requiring specific data declare trait bounds. Missing data is a compile error.

### Traversal

Traversal methods are built into `Graph`:

- `par_iter_breadth_first_forward` / `par_iter_breadth_first_backward` ([packages/treetime-graph/src/graph_traverse.rs#L221-L250](../../packages/treetime-graph/src/graph_traverse.rs#L221-L250)) - parallel breadth-first using rayon, processing independent nodes concurrently within each frontier
- `iter_depth_first_preorder_forward` / `iter_depth_first_postorder_forward` ([packages/treetime-graph/src/graph_traverse.rs#L256-L321](../../packages/treetime-graph/src/graph_traverse.rs#L256-L321)) - sequential depth-first
- `iter_breadth_first_forward` / `iter_breadth_first_reverse` ([packages/treetime-graph/src/graph_traverse.rs#L327-L354](../../packages/treetime-graph/src/graph_traverse.rs#L327-L354)) - sequential breadth-first

Traversal callbacks receive `GraphNodeForward` (`#GraphNodeForward`) or `GraphNodeBackward` (`#GraphNodeBackward`) in [packages/treetime-graph/src/graph_traverse.rs#L18-L167](../../packages/treetime-graph/src/graph_traverse.rs#L18-L167), providing mutable payload access together with read access to parent/child payloads. `GraphNodeForward` gives mutable access to child edges with read-only parents (natural for root-to-leaf passes). `GraphNodeBackward` gives mutable access to parent edges with read-only children (natural for leaf-to-root passes).

### Parallel BFS engine

The parallel BFS engine in [packages/treetime-graph/src/breadth_first.rs#L139-L201](../../packages/treetime-graph/src/breadth_first.rs#L139-L201) uses level-synchronous traversal with rayon:

1. Start with source nodes (roots or leaves) as the initial frontier.
2. Process each frontier node in parallel via `rayon::into_par_iter()`.
3. For each unvisited node: call the explorer callback, mark as visited, collect successors.
4. Filter successors: add to next frontier only if ALL predecessors are already visited.
5. Repeat until the frontier is empty.

The predecessor-checking logic at step 4 handles nodes with multiple inbound edges (multiple parents in a network). A node waits until all its predecessors are processed before entering the frontier. This is the standard level-synchronous BFS pattern (Buluc & Madduri, 2011 [[29]](#ref-29)) applied to DAGs, where the number of synchronization barriers equals the longest path length in the graph. For balanced trees this is O(log n). The parallel traversal methods start from `self.roots` (plural) without enforcing a single root.

### Topology manipulation

The graph provides primitive topology operations in [packages/treetime-graph/src/graph_ops.rs](../../packages/treetime-graph/src/graph_ops.rs):

- `add_node(payload) -> GraphNodeKey` - append a node, return its key
- `add_edge(source, target, payload) -> GraphEdgeKey` - create a directed edge; validates no self-loops and no duplicate parallel edges
- `remove_node(key)` - remove a node and all incident edges
- `remove_edge(key)` - remove an edge, update both endpoints' adjacency lists
- `collapse_edge(edge_key)` - merge the target node into the source node, redirecting all of the target's other edges to the source
- `build()` - recompute cached root and leaf sets

The `find_paths` module in [packages/treetime-graph/src/find_paths.rs](../../packages/treetime-graph/src/find_paths.rs) provides DAG-aware path queries: `find_paths()` returns all edges on all paths between two nodes, and `exists_forward_path_between()` / `exists_backward_path_between()` check path existence via BFS with early termination. These operate on general directed graphs without assuming tree topology, making them building blocks for network algorithms such as finding displayed trees or identifying reticulation cycles.

The `invert_edge()` (`#invert_edge`) free function in [packages/treetime-graph/src/edge.rs#L102-L133](../../packages/treetime-graph/src/edge.rs#L102-L133) reverses an edge's direction by swapping source/target and updating both nodes' inbound/outbound lists. This is used during rerooting: `apply_reroot()` (`#apply_reroot`) in [packages/treetime/src/commands/clock/reroot.rs#L196-L228](../../packages/treetime/src/commands/clock/reroot.rs#L196-L228) finds the path from the new root to the old root and inverts every edge along that path.

### Tree enforcement

The storage layer is a general directed graph. The `add_edge` validation prevents self-loops and duplicate parallel edges but does not prevent multiple edges pointing to the same target from different sources (multiple parents). Multiple roots are representable.

Tree semantics are enforced at runtime by methods that return errors or panic when the graph is not a tree:

- `one_parent_of()` (`#one_parent_of`) in [packages/treetime-graph/src/graph.rs#L121-L137](../../packages/treetime-graph/src/graph.rs#L121-L137) returns `Ok(None)` for roots, `Ok(Some(...))` for single parent, and an error for multiple parents.
- `get_exactly_one_root()` in [packages/treetime-graph/src/graph.rs#L261-L271](../../packages/treetime-graph/src/graph.rs#L261-L271) returns an error if there is not exactly one root.
- `GraphNodeForward::get_exactly_one_parent()` in [packages/treetime-graph/src/graph_traverse.rs#L86-L90](../../packages/treetime-graph/src/graph_traverse.rs#L86-L90) errors on multiple parents.
- `GraphNodeBackward::get_exactly_one_parent_edge()` in [packages/treetime-graph/src/graph_traverse.rs#L165-L167](../../packages/treetime-graph/src/graph_traverse.rs#L165-L167) errors on multiple parents.

All sequential traversal methods call `get_exactly_one_root().unwrap()` and panic on multi-root graphs. All algorithm code in v1 currently calls one of these enforcement methods, making v1 strictly tree-only at the algorithm level.

The test `test_collapse_edge_multiple_inbound_edges` in [packages/treetime/src/graph/\_\_tests\_\_/graph.rs#L303-L351](../../packages/treetime/src/graph/__tests__/graph.rs#L303-L351) explicitly constructs a graph where one node has two parents and verifies that `collapse_edge()` handles it correctly. This confirms the storage and mutation layer supports network topology even though algorithms do not use it.

## Rationale

### Immediate requirements

**No typed payloads.** Algorithm data in v0 is stored as dynamic attributes on Clade objects. Any attribute can be set on any node at any time. Missing attributes cause runtime `AttributeError`. v1 uses Rust's type system to ensure each algorithm stage has exactly the data it needs at compile time.

**No explicit edges.** In Bio.Phylo, branch data (length, mutations) is stored on the child node (BioPython docs [[28]](#ref-28)). There is no edge entity. v1's `Edge<E>` provides a natural place for branch-specific data like length distributions and message-passing fields.

**No parent pointers.** BioPython trees are child-linked only (BioPython docs [[28]](#ref-28)). v0 patches parent references manually in `_prepare_nodes()`. v1 stores bidirectional edge references (inbound/outbound) as part of the node structure, making parent access a built-in O(1) operation via `node.inbound()` (`#inbound`) in [packages/treetime-graph/src/node.rs#L190-L192](../../packages/treetime-graph/src/node.rs#L190-L192).

**Single-threaded traversal.** BioPython's `find_clades` yields nodes sequentially. v1's parallel breadth-first traversal processes independent nodes concurrently, using frontier-based synchronization to respect data dependencies between tree levels.

### Forward-looking design for phylogenetic networks

The directed graph representation is deliberately more general than the tree algorithms currently require. Several design choices prepare for future network support:

**Bidirectional edge storage.** Each node stores both `inbound_edges` and `outbound_edges`, allowing any node to have multiple parents. A child-linked tree representation would require structural changes to support reticulate nodes. The directed graph representation requires only relaxing runtime checks.

**Edge as first-class entity.** In phylogenetic networks, the edge connecting a reticulate node to each parent carries distinct data: branch length, inheritance probability, and its own substitution model. Storing branch data on the child node (as in BioPython) cannot distinguish between contributions from different parents. v1's explicit edge entity handles this naturally.

**Multiple root support.** The `roots` vector holds multiple root keys. ARGs and some phylogenetic network representations have multiple roots (one per segment or per region of the genome). The storage layer and parallel BFS already iterate over all roots. Only the sequential traversal methods and algorithm code enforce single root.

**DAG-aware traversal.** The parallel BFS engine checks that all predecessors of a node are visited before processing it. This correctly handles reticulate nodes where two parent lineages must be processed before the child. A tree-only BFS would not need this check (each node has exactly one predecessor).

**Softwired parsimony path.** Extending Fitch's parsimony to networks requires selecting one parent edge per reticulate node per site. The 2-approximation algorithm for softwired parsimony on tree-child networks (Frohn & Kelk, 2024 [[16]](#ref-16)) uses a modified Fitch pass where reticulate nodes inherit state from the selected parent. v1's traversal context types (`GraphNodeForward`, `GraphNodeBackward`) already expose all parent/child payloads, providing the data access pattern needed for this extension.

**Level-synchronous parallelism on networks.** For level-k networks with bounded treewidth (Markin et al., 2024 [[12]](#ref-12)), the BFS frontier approach remains efficient: each frontier processes nodes at the same DAG depth, and the number of barriers equals the longest directed path. The frontier size may increase at reticulate nodes (a node appears in the frontier only when all its predecessors are done), but the algorithm is unchanged.

## Practical impact

- Node data access is type-checked at compile time. Algorithms that require specific fields (e.g., time distributions) declare trait bounds, and missing data is a compile error rather than a runtime `AttributeError`.
- Branch data lives on edges. Functions that operate on branches receive `Edge<E>` payloads directly instead of reading attributes from child nodes.
- Traversal order is explicit in the method name. There is no string parameter (`order='postorder'`) that can be misspelled.
- Parallel traversal reduces wall-clock time for large trees. The BFS frontier approach processes all nodes at the same tree depth concurrently.
- The graph structure accepts phylogenetic networks. Trees are a special case where each node has at most one parent. Algorithms enforce tree semantics at runtime via `one_parent_of()` and `get_exactly_one_root()`, localizing the constraint to specific methods that can be relaxed as network algorithms are implemented.

## Relevance to viral phylodynamics

TreeTime's primary application is viral phylodynamics: inferring evolutionary rates, divergence times, and transmission dynamics from pathogen genome sequences (Sagulenko et al., 2018 [[34]](#ref-34)). Viral evolution is precisely the domain where the tree assumption breaks down most frequently.

SARS-CoV-2 recombinant lineages (XBB, XEC, XFG) arise when a host is co-infected with two distinct variants and template switching produces a chimeric genome (Wikipedia [[7]](#ref-7)). A TreeTime analysis of a recombinant lineage assigns it to a single branch in the tree, but different genomic regions of the same virus have different phylogenetic histories. The resulting tree is a compromise that does not reflect the true evolutionary path of any region. Network-aware inference would represent the recombination event as a reticulate node with two parent edges, one for each parental region.

Influenza reassortment creates the same problem at the segment level. When TreeTime builds a timetree from concatenated influenza segments, a reassortant virus appears on a single branch. A network representation would place the reassortant at a hybrid node whose parent edges connect to the donor lineages for each segment, with the breakpoint at segment boundaries rather than within a continuous genome.

These scenarios are not edge cases. Recombination detection is a routine step in SARS-CoV-2 genomic surveillance, and reassortment tracking is central to influenza pandemic preparedness. A graph representation that can encode reticulate nodes positions TreeTime to move from detecting these events as anomalies to modeling them as part of the evolutionary history.

## Related software

**SplitsTree** (Huson & Bryant, 2006 [[30]](#ref-30)). Java, University of Tubingen. The primary tool for unrooted network visualization. Implements split decomposition, NeighborNet, consensus networks, super networks, and median networks. Input via NEXUS files with a SPLITS data block. Produces planar split graphs where incompatible splits appear as parallelogram-shaped regions. https://splitstree.org

**PhyloNet** (Wen et al., 2018 [[37]](#ref-37)). Java, Rice University (Nakhleh Lab). The most complete toolkit for rooted network inference under the multispecies network coalescent. Inference methods include maximum parsimony (`InferNetwork_MP`), maximum likelihood (`InferNetwork_ML`), maximum pseudo-likelihood (`InferNetwork_MPL`), and Bayesian MCMC from gene trees (`MCMC_GT`) or directly from sequence alignments (`MCMC_SEQ`). Input/output uses Rich Newick format (Extended Newick with inheritance probabilities). Supports simulation (`SimGTinNetwork`), network comparison (`Cmpnets`), and introgression quantification (`CallIntroRate`). https://phylogenomics.rice.edu/html/phylonet.html

**PhyloNetworks** (Solis-Lemus & Ane, 2016 [[19]](#ref-19)). Julia, University of Wisconsin. Implements SNaQ (Species Networks applying Quartets), a pseudolikelihood approach using quartet concordance factors under the coalescent model with hybridization. Scales to hundreds of taxa. Also supports trait evolution (Brownian motion, Pagel's lambda) on networks via belief propagation. https://github.com/crsl4/PhyloNetworks.jl

**Dendroscope** (Huson & Bryant, 2006 [[30]](#ref-30)). Java, University of Tubingen. Interactive viewer for phylogenetic trees and networks of all sizes. Renders Extended Newick networks as DAGs. Does not display inheritance probabilities from Rich Newick. https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/dendroscope/

**SpeciesNetwork/BEAST 2** (Zhang et al., 2018 [[20]](#ref-20)). Java, BEAST 2 add-on. Bayesian inference of time-calibrated species networks from multilocus sequence data. Uses a birth-hybridization process prior for the network and the MSNC prior for embedded gene trees. Supports relaxed molecular clock models. Jointly estimates network topology, divergence times, population sizes, and inheritance probabilities via MCMC.

**tskit** (Kelleher et al., 2018 [[35]](#ref-35)). Python/C, University of Oxford. Efficient data structure and library for storing and manipulating tree sequences (succinct representations of correlated local trees along a genome). Provides the foundation for population-scale ARG inference and simulation. Used by msprime (coalescent simulator), SLiM (forward simulator), and tsinfer/tsdate (ARG inference from sequence data). https://tskit.dev

## Books

- Huson, D.H., Rupp, R. & Scornavacca, C. (2011). **Phylogenetic Networks: Concepts, Algorithms and Applications** [[4]](#ref-4). Cambridge University Press. The primary textbook on the subject. Covers split networks, median networks, Buneman graphs, NeighborNet, consensus and super networks, hybridization networks, ARGs, rooted network classes (galled trees, level-k, tree-child, tree-sibling, tree-based), tight span theory, and distance matrix methods. Each chapter pairs mathematical foundations with algorithmic implementations.
- Semple, C. & Steel, M. (2003). **Phylogenetics** [[38]](#ref-38). Oxford University Press. Mathematical foundations of phylogenetic trees: splits, compatibility, tree metrics, reconstruction algorithms, statistical consistency. Prerequisite for the network extensions in Huson et al.
- Steel, M. (2016). **Phylogeny: Discrete and Random Processes in Evolution** [[39]](#ref-39). SIAM. Probabilistic and combinatorial theory of phylogenetic trees and networks, including Markov models on trees, identifiability, and ancestral reconstruction consistency results.
- Gontier, N. (2015). **Reticulate Evolution** [[5]](#ref-5). Springer. Biological foundations: hybridization, HGT, symbiogenesis, and their phylogenetic implications across all domains of life.

## Further topics

The following topics are tangential to the v0/v1 deviation but provide depth for readers extending v1 toward network support.

### Inheritance probabilities

At a reticulate node with two parents, the inheritance probability gamma (0 < gamma < 1) specifies the fraction of the genome (or the probability for any given locus) inherited from one parent, with 1 - gamma from the other. Under the multispecies network coalescent, a gene lineage entering a reticulate node follows parent edge A with probability gamma and parent edge B with probability 1 - gamma, then coalesces independently within each population (Elworth et al., 2019 [[18]](#ref-18)). PhyloNet encodes gamma as the third colon-delimited field in Rich Newick: `(child)parent#H1:branch_length:population_size:gamma`. Inheritance probabilities are a key parameter in likelihood-based network inference and must be stored on edges, not nodes, since each parent edge carries a different gamma.

### Endosymbiosis

Endosymbiosis is an extreme form of reticulation where one organism lives inside another, leading to permanent genetic merger. The mitochondrion originated from an alpha-proteobacterial endosymbiont approximately 2 billion years ago; the chloroplast from a cyanobacterial endosymbiont approximately 1.5 billion years ago (Gontier, 2015 [[5]](#ref-5)). Over time, genes transferred from the endosymbiont genome to the host nuclear genome (endosymbiotic gene transfer). The eukaryotic tree of life is a network at its root: the host and endosymbiont lineages merge into the ancestral eukaryote. Secondary and tertiary endosymbioses (an alga engulfed by another eukaryote) created further reticulation in protist and algal phylogenies.

### Tree containment

The tree containment problem asks: given a phylogenetic network N and a phylogenetic tree T, does N display T? That is, can T be obtained from N by selecting one parent edge per reticulate node and suppressing degree-2 nodes? This problem is polynomial for tree-child, tree-sibling, and level-k networks (for fixed k), but NP-complete for general phylogenetic networks (Huson et al., 2011 [[4]](#ref-4)). Tree containment is relevant to validation: given an inferred network and a set of gene trees, checking whether each gene tree is displayed by the network tests consistency of the network with the data.

### Buneman graphs and median networks

The Buneman graph of a set of taxa and characters places each taxon at a vertex and adds a latent vertex for every combination of character values where every pairwise combination is observed in the data. Every maximum parsimony tree for the character data embeds as a subgraph of the Buneman graph. The Buneman graph is a tree if and only if a perfect phylogeny exists (no pair of incompatible characters with all four value combinations observed) (Huson et al., 2011 [[4]](#ref-4)). Median networks are a special case for binary characters: given a set of binary sequences, the median network contains all Steiner trees (minimum spanning networks) connecting the input sequences. Median networks are widely used in human mitochondrial DNA and Y-chromosome phylogeography, where low mutation rates make tree-like representations appropriate but recurrent mutations create localized incompatibilities.

### Traversal fusion

When multiple traversals over the same tree are needed (Fitch parsimony pass, marginal backward pass, marginal forward pass, time optimization pass), fusing them into fewer passes reduces total work. Sakka et al. (2019 [[40]](#ref-40)) describe fine-grained traversal fusion for heterogeneous trees, where different node types (leaf vs internal) carry different data. The technique merges compatible traversals at the individual node level rather than requiring whole passes to have identical structure. For phylogenetic trees where leaf nodes and internal nodes have different payload types, this permits fusing an upward likelihood pass with an upward parsimony pass even though they operate on different data. v1's typed payload system (`NodeAncestral` vs `NodeTimetree`) maps to the heterogeneous tree model. Traversal fusion is orthogonal to network support but compounds with it: network algorithms require more passes than tree algorithms, making fusion more valuable.

## References

<a id="ref-1"></a>[1] Qi, J. & Schicho, J. (2020). Five equivalent representations of a phylogenetic tree. arXiv:2011.11774.

<a id="ref-2"></a>[2] Felsenstein, J. The Newick tree format. PHYLIP documentation. https://phylipweb.github.io/phylip/newicktree.html

<a id="ref-3"></a>[3] Talevich, E., Invergo, B.M., Cock, P.J. & Chapman, B.A. (2012). Bio.Phylo: A unified toolkit for processing, analyzing and visualizing phylogenetic trees in Biopython. BMC Bioinformatics 13:209.

<a id="ref-4"></a>[4] Huson, D.H., Rupp, R. & Scornavacca, C. (2011). Phylogenetic Networks: Concepts, Algorithms and Applications. Cambridge University Press.

<a id="ref-5"></a>[5] Gontier, N. (2015). Reticulate Evolution. Springer.

<a id="ref-6"></a>[6] Ge, F., Wang, L.S. & Kim, J. (2005). The cobweb of life revealed by genome-scale estimates of horizontal gene transfer. PLOS Biology 3:e316.

<a id="ref-7"></a>[7] Wikipedia contributors. SARS-CoV-2. https://en.wikipedia.org/wiki/SARS-CoV-2

<a id="ref-8"></a>[8] Wikipedia contributors. Reassortment. https://en.wikipedia.org/wiki/Reassortment

<a id="ref-9"></a>[9] Wikipedia contributors. Incomplete lineage sorting. https://en.wikipedia.org/wiki/Incomplete_lineage_sorting

<a id="ref-10"></a>[10] Huson, D.H. & Scornavacca, C. (2011). A survey of combinatorial methods for phylogenetic networks. Genome Biology and Evolution 3:23-35.

<a id="ref-11"></a>[11] Pons, J.C., Semple, C. & Steel, M. (2018). Tree-based networks: characterisations, metrics, and support trees. Journal of Mathematical Biology 78:899-918. arXiv:1710.07836.

<a id="ref-12"></a>[12] Markin, A. et al. (2024). Treewidth bounds for level-k phylogenetic networks. arXiv:2411.13380.

<a id="ref-13"></a>[13] Fischer, M., van Iersel, L., Kelk, S. & Scornavacca, C. (2015). On Computing the Maximum Parsimony Score of a Phylogenetic Network. SIAM Journal on Discrete Mathematics 29:559-585. arXiv:1302.2430.

<a id="ref-14"></a>[14] Wikipedia contributors. Split (phylogenetics). https://en.wikipedia.org/wiki/Split_(phylogenetics)

<a id="ref-15"></a>[15] Wikipedia contributors. Coalescent theory. https://en.wikipedia.org/wiki/Coalescent_theory

<a id="ref-16"></a>[16] Frohn, M. & Kelk, S. (2024). A 2-approximation algorithm for the softwired parsimony problem on binary, tree-child phylogenetic networks. arXiv:2409.18077.

<a id="ref-17"></a>[17] Doecker, C., Linz, S. & Wicke, K. (2024). Bounding the softwired parsimony score of a phylogenetic network. arXiv:2405.19587.

<a id="ref-18"></a>[18] Elworth, R.A.L. et al. (2019). Advances in Computational Methods for Phylogenetic Networks in the Presence of Hybridization. arXiv:1808.08662.

<a id="ref-19"></a>[19] Solis-Lemus, C. & Ane, C. (2016). Inferring Phylogenetic Networks with Maximum Pseudolikelihood under Incomplete Lineage Sorting. PLOS Genetics 12:e1005896.

<a id="ref-20"></a>[20] Zhang, C. et al. (2018). Bayesian Inference of Species Networks from Multilocus Sequence Data. Molecular Biology and Evolution 35:504-517.

<a id="ref-21"></a>[21] Xu, W. & Ane, C. (2024). A consistent least-squares criterion for calibrating edge lengths in phylogenetic networks. arXiv:2407.19343.

<a id="ref-22"></a>[22] Mitchell, J.D. & Holland, B.R. (2025). Convergence-divergence models: Generalizations of phylogenetic trees modeling gene flow over time. arXiv:2504.07384.

<a id="ref-23"></a>[23] Francis, A. & Moulton, V. (2018). Identifiability of tree-child phylogenetic networks under a probabilistic recombination-mutation model of evolution. Journal of Theoretical Biology 446:160-167. arXiv:1712.04223.

<a id="ref-24"></a>[24] Allman, E.S., Ane, C., Banos, H. & Rhodes, J.A. (2025). Identifiability of galled tree-child networks from quartet concordance factors. arXiv:2504.21116.

<a id="ref-25"></a>[25] van Iersel, L. & Moulton, V. (2014). Trinets encode tree-child and level-2 phylogenetic networks. Journal of Mathematical Biology 68:1707-1729. arXiv:1210.0362.

<a id="ref-26"></a>[26] To, T.H. & Habib, M. (2009). Level-k phylogenetic networks are constructable from a dense triplet set in polynomial time. arXiv:0901.1657.

<a id="ref-27"></a>[27] Cardona, G., Rossello, F. & Valiente, G. (2008). Extended Newick: it is time for a standard representation of phylogenetic networks. BMC Bioinformatics 9:532. PMC2621367.

<a id="ref-28"></a>[28] BioPython Bio.Phylo.BaseTree API documentation. https://biopython.org/docs/latest/api/Bio.Phylo.BaseTree.html

<a id="ref-29"></a>[29] Buluc, A. & Madduri, K. (2011). Parallel Breadth-First Search on Distributed Memory Systems. arXiv:1104.4518.

<a id="ref-30"></a>[30] Huson, D.H. & Bryant, D. (2006). Application of phylogenetic networks in evolutionary studies. Molecular Biology and Evolution 23:254-267.

<a id="ref-31"></a>[31] Green, R.E. et al. (2010). A draft sequence of the Neandertal genome. Science 328:710-722.

<a id="ref-32"></a>[32] Wikipedia contributors. Gene duplication. https://en.wikipedia.org/wiki/Gene_duplication

<a id="ref-33"></a>[33] Bryant, D. & Moulton, V. (2004). Neighbor-Net: An Agglomerative Method for the Construction of Phylogenetic Networks. Molecular Biology and Evolution 21:255-265.

<a id="ref-34"></a>[34] Sagulenko, P., Puller, V. & Neher, R.A. (2018). TreeTime: Maximum-likelihood phylodynamic analysis. Virus Evolution 4:vex042.

<a id="ref-35"></a>[35] Kelleher, J., Thornton, K.R., Ashander, J. & Ralph, P.L. (2018). Efficient pedigree recording for fast population genetics simulation. PLOS Computational Biology 14:e1006581. https://doi.org/10.1371/journal.pcbi.1006581

<a id="ref-36"></a>[36] Fournier, R. & Larribe, F. (2025). Moonshine.jl: a Julia package for genome-scale model-based ancestral recombination graph inference. arXiv:2511.21124.

<a id="ref-37"></a>[37] Wen, D., Yu, Y., Zhu, J. & Nakhleh, L. (2018). Inferring phylogenetic networks using PhyloNet. Systematic Biology 67:735-740.

<a id="ref-38"></a>[38] Semple, C. & Steel, M. (2003). Phylogenetics. Oxford University Press.

<a id="ref-39"></a>[39] Steel, M. (2016). Phylogeny: Discrete and Random Processes in Evolution. SIAM.

<a id="ref-40"></a>[40] Sakka, L., Krishnamurthy, R., Kuber, K., Ragan-Kelley, J. & Newton, R.R. (2019). Sound, Fine-Grained Traversal Fusion for Heterogeneous Trees. arXiv:1904.07061.
