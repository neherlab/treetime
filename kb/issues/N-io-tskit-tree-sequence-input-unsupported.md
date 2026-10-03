# tskit tree sequence input is not supported

TreeTime cannot read tskit tree sequences (`.trees`). sc2ts stores an ancestral recombination graph of 2.48 million SARS-CoV-2 genomes in this format ([Zhan et al. 2023](https://doi.org/10.1101/2023.06.08.544212)). The TreeTime graph already supports networks ([kb/decisions/graph-based-phylogenetic-representation.md](../decisions/graph-based-phylogenetic-representation.md)), so reading the graph structure is possible in principle.

## Open questions

- Scope: does TreeTime need time inference on recombinant networks, and do the inference algorithms support nodes with more than one parent?
- A tree sequence assigns each edge to a genome interval. How do these intervals map onto partitions?
- Library: the tech stack lists no tskit library, so reading the format needs a dependency decision

## Related

- [kb/proposals/io-format-coverage.md](../proposals/io-format-coverage.md): ranking of format work
