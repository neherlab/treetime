# tskit tree sequence input is not supported

TreeTime cannot read tskit tree sequences (`.trees`). sc2ts stores an ancestral recombination graph of 2.48 million SARS-CoV-2 genomes in this format ([Zhan et al. 2023](https://doi.org/10.1101/2023.06.08.544212)). The TreeTime graph already supports networks ([kb/decisions/graph-based-phylogenetic-representation.md](../decisions/graph-based-phylogenetic-representation.md)), so reading the graph structure is possible in principle.

The Newick reader rejects every node with more than one parent: `fn graph_from_newick()` calls `fn reject_network_nodes()`, which stops with "TreeTime reads trees only" [packages/treetime-io/src/nwk.rs#L48-L49](../../packages/treetime-io/src/nwk.rs#L48-L49), [packages/treetime-io/src/nwk.rs#L105-L118](../../packages/treetime-io/src/nwk.rs#L105-L118).

## Open questions

> [!IMPORTANT]
> **Decision required.** The reader needs three choices before work can start:
>
> - Scope: does TreeTime need time inference on recombinant networks? Options: extend the inference algorithms to nodes with more than one parent, or read a tree sequence as a set of local trees
> - Partitions: a tree sequence assigns each edge to a genome interval. It is open how these intervals map onto partitions
> - Library: the tech stack lists no tskit library. Options: add a maintained Rust tskit library through the dependency process, or read the format with a project reader

## Related

- [kb/proposals/io-format-coverage.md](../proposals/io-format-coverage.md): ranking of format work
