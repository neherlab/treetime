# Graph-based phylogenetic representation

TreeTime v1 represents phylogenetic structure as a directed graph with explicit nodes and edges. A phylogenetic tree is the restricted case in which every non-root node has exactly one parent. The representation can also express multiple-parent ancestry, while current inference algorithms continue to require tree topology.

## Requirements from biology

Evolutionary histories include several relationship patterns:

- **Divergence:** one ancestral lineage gives rise to descendant lineages. A rooted tree represents this directly.
- **<a id="gloss-use-1"></a>Reticulation <sup>[1](#gloss-1)</sup>:** hybridization, introgression, horizontal gene transfer, viral recombination, and reassortment combine ancestry from multiple lineages. The resulting history requires nodes with multiple parents.
- **Locus-dependent ancestry:** recombination and reassortment can give different genomic regions different histories. A representation may need to associate an ancestry path with a genomic interval or segment.
- **Conflicting signal:** incompatible characters or distances may show that no single tree fits the observations, even when the data do not identify a particular reticulate history.
- **Edge-specific information:** branch length, mutations, inheritance contribution, and genomic scope describe a relationship between two nodes. They belong to edges rather than either endpoint alone.

Phylogenetic-network methods distinguish explicit evolutionary histories from visualizations of incompatible signal, and impose different structural restrictions for different biological and computational settings <a id="cite-1"></a>[Huson and Scornavacca 2011](https://doi.org/10.1093/gbe/evq077) [[1](#ref-1)]. TreeTime therefore needs a representation that preserves directed ancestry without treating every network formalism as equivalent.

## Representation alternatives

| Representation | What it expresses | Fit for TreeTime |
| :--- | :--- | :--- |
| Rooted tree | One directed ancestry path from the root to each sample | Matches current inference algorithms, but cannot represent reticulation |
| Collection of local or <a id="gloss-use-2"></a>displayed trees <sup>[2](#gloss-2)</sup> | A different tree for each locus, segment, or parent-edge selection | Represents locus-specific outcomes, but does not preserve reticulation events as shared entities |
| <a id="gloss-use-3"></a>Split network <sup>[3](#gloss-3)</sup> | Incompatible splits or distances in an unrooted graph | Useful for exploratory signal visualization, but does not necessarily describe directed evolutionary history |
| Rooted phylogenetic network | Divergence and reticulation in a rooted <a id="gloss-use-4"></a>directed acyclic graph (DAG) <sup>[4](#gloss-4)</sup> | Directly generalizes a rooted tree and preserves explicit ancestry paths |
| <a id="gloss-use-5"></a>Ancestral recombination graph (ARG) <sup>[5](#gloss-5)</sup> or <a id="gloss-use-6"></a>tree sequence <sup>[6](#gloss-6)</sup> | Changing genealogies along a genome, including the intervals inherited through each edge | Appropriate when genomic coordinates are intrinsic to the model; requires semantics beyond generic topology |

ARG work makes the additional requirement explicit: adjacent nucleotides can follow different inheritance paths, so a complete recombination representation associates genomes and genomic intervals with ancestry <a id="cite-2"></a>[Wong et al. 2024](https://doi.org/10.1093/genetics/iyae100) [[2](#ref-2)]. A generic directed graph can host those annotations, but does not define their meaning by itself.

## Decision

Use a rooted directed graph as TreeTime's common phylogenetic representation:

- Nodes represent observed samples or inferred ancestors.
- Directed edges represent ancestry from parent to child and own branch-specific data.
- Trees remain a first-class restricted topology rather than a separate data model.
- Multiple parents and multiple roots are representable when an input format and algorithm define their semantics.
- Node, edge, and graph payloads are typed so algorithms declare the biological data they require.

This choice captures the common topology shared by trees and explicit rooted phylogenetic networks. It also lets tree algorithms state their one-parent and one-root preconditions at the algorithm boundary.

## Why this representation

### Direct correspondence with ancestry

Directed edges preserve the order of descent. A reticulate node can retain each contributing parent relationship explicitly, unlike a tree node that can name only one parent. Forward and backward traversals can respect ancestry dependencies without making traversal bookkeeping part of the biological model.

### Explicit branch identity

An edge is a biological entity, not an attribute of its child node. This permits each parent contribution to carry distinct branch length, mutations, inheritance probability, substitution parameters, or genomic scope.

### One model for current and broader topology

TreeTime's existing operations use tree topology. Representing trees as constrained graphs avoids conversion between unrelated tree and network object models, while keeping topology requirements visible in each algorithm's contract.

### Extensible semantics

The graph supplies identity, direction, and connectivity. Specialized models can add interval annotations for ARGs, inheritance probabilities for reticulate edges, or segment identifiers for reassortment without changing the core topological concept.

## Boundaries of the decision

The graph representation does not by itself provide:

- proof that a topology is acyclic or biologically valid;
- a likelihood, parsimony, molecular-clock, or inheritance model for networks;
- genomic interval semantics required by ARGs and tree sequences;
- an interpretation of incompatible signal as a specific evolutionary event;
- serialization of reticulate identities or inheritance annotations.

These belong to validated input types, algorithm contracts, and format adapters. In particular, split networks and rooted phylogenetic networks answer different questions: one can visualize incompatible evidence, while the other proposes an explicit directed history.

## Ecosystem precedents

Other phylogenetic systems choose representations according to their scientific scope:

- **SplitsTree** combines tree, implicit-network, and explicit-network methods, illustrating that conflict visualization and explicit reticulation are distinct representations <a id="cite-3"></a>[Kloepper and Huson 2008](https://doi.org/10.1186/1471-2148-8-22) [[3](#ref-3)].
- **PhyloNet** represents and infers explicit phylogenetic networks from multilocus or sequence data using parsimony, likelihood, Bayesian, and pseudolikelihood methods <a id="cite-4"></a>[Wen et al. 2018](https://doi.org/10.1093/sysbio/syy015) [[4](#ref-4)].
- **PhyloNetworks** uses explicit networks for pseudolikelihood inference that accounts for incomplete lineage sorting and reticulation <a id="cite-5"></a>[Solís-Lemus and Ané 2016](https://doi.org/10.1371/journal.pgen.1005896) [[5](#ref-5)].
- **tskit** uses succinct tree sequences to retain correlated local genealogies and their genomic scope efficiently <a id="cite-6"></a>[Kelleher et al. 2018](https://doi.org/10.1371/journal.pcbi.1006581) [[6](#ref-6)].

TreeTime focuses on rapid phylodynamic analysis of rooted phylogenies <a id="cite-7"></a>[Sagulenko et al. 2018](https://doi.org/10.1093/ve/vex042) [[7](#ref-7)]. A directed graph is the appropriate common representation because it preserves that rooted workflow and admits explicit reticulation, without claiming the additional semantics of split networks, ARGs, or network-inference models.

## Current support

- TreeTime's inference algorithms enforce tree topology.
- The graph model can represent nodes with multiple parents and more than one root.
- Traversal supports both ancestry directions and processes a node only after its dependencies in that direction are complete.
- Standard Newick input and output remain tree-oriented; network serialization requires a format that preserves reticulate node identity and edge annotations.
- Algorithm and source-code details are maintained in [kb/algo/graph.md](../algo/graph.md) and the `treetime-graph` crate rather than duplicated in this architectural decision.

## Consequences

- Tree algorithms retain explicit and testable topology preconditions.
- Branch-specific biological data has an unambiguous owner.
- Network-aware algorithms can reuse graph identity and traversal concepts while defining their own scientific semantics.
- Importers must distinguish a directed graph that is structurally representable from a phylogenetic network that is scientifically valid for a particular analysis.
- Supporting ARGs or tree sequences requires interval-aware data and operations in addition to graph topology.

## Glossary

1. <a id="gloss-1"></a> **Reticulation.** An evolutionary event in which a lineage receives ancestry from more than one parent lineage. [↩](#gloss-use-1)
2. <a id="gloss-2"></a> **Displayed tree.** A tree obtained from a rooted phylogenetic network by selecting one incoming edge at each reticulate node and suppressing redundant intermediate nodes. [↩](#gloss-use-2)
3. <a id="gloss-3"></a> **Split network.** An unrooted network that represents compatible or incompatible bipartitions of a taxon set, commonly for visualizing conflicting signal. [↩](#gloss-use-3)
4. <a id="gloss-4"></a> **Directed acyclic graph (DAG).** A directed graph containing no directed cycle. Rooted phylogenetic networks use this structure to preserve ancestry order. [↩](#gloss-use-4)
5. <a id="gloss-5"></a> **Ancestral recombination graph (ARG).** A representation of genome ancestry in which recombination permits different genomic intervals to follow different genealogical paths. [↩](#gloss-use-5)
6. <a id="gloss-6"></a> **Tree sequence.** A compact representation of correlated local genealogical trees and the genomic intervals over which their edges apply. [↩](#gloss-use-6)

## References

1. <a id="ref-1"></a> Huson, Daniel H., and Celine Scornavacca. 2011. "A survey of combinatorial methods for phylogenetic networks." _Genome Biology and Evolution_ 3:23-35. https://doi.org/10.1093/gbe/evq077 [↩](#cite-1)
2. <a id="ref-2"></a> Wong, Yan, Anastasia Ignatieva, Jere Koskela, Gregor Gorjanc, Anthony W. Wohns, and Jerome Kelleher. 2024. "A general and efficient representation of ancestral recombination graphs." _Genetics_ 228:iyae100. https://doi.org/10.1093/genetics/iyae100 [↩](#cite-2)
3. <a id="ref-3"></a> Kloepper, Tobias H., and Daniel H. Huson. 2008. "Drawing explicit phylogenetic networks and their integration into SplitsTree." _BMC Evolutionary Biology_ 8:22. https://doi.org/10.1186/1471-2148-8-22 [↩](#cite-3)
4. <a id="ref-4"></a> Wen, Dingqiao, Yun Yu, Jiafan Zhu, and Luay Nakhleh. 2018. "Inferring phylogenetic networks using PhyloNet." _Systematic Biology_ 67:735-740. https://doi.org/10.1093/sysbio/syy015 [↩](#cite-4)
5. <a id="ref-5"></a> Solís-Lemus, Claudia, and Cécile Ané. 2016. "Inferring phylogenetic networks with maximum pseudolikelihood under incomplete lineage sorting." _PLOS Genetics_ 12:e1005896. https://doi.org/10.1371/journal.pgen.1005896 [↩](#cite-5)
6. <a id="ref-6"></a> Kelleher, Jerome, Kevin R. Thornton, Jaime Ashander, and Peter L. Ralph. 2018. "Efficient pedigree recording for fast population genetics simulation." _PLOS Computational Biology_ 14:e1006581. https://doi.org/10.1371/journal.pcbi.1006581 [↩](#cite-6)
7. <a id="ref-7"></a> Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-likelihood phylodynamic analysis." _Virus Evolution_ 4:vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-7)
