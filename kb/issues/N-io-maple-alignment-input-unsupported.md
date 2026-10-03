# MAPLE alignment input is not supported

TreeTime cannot read alignments in the MAPLE format. The format stores a reference sequence once and then the differences of each sample from that reference, so it is much smaller than FASTA for large sets of closely related genomes. CMAPLE reads it [[src](https://github.com/iqtree/cmaple/blob/3d45b1ab68e2d68a2825bf17a531e22200578cd6/alignment/alignment.h#L15-L18)], and IQ-TREE 3 can write it [[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/alignment/alignment.cpp#L6297)].

The per-sample difference lists map directly onto sparse partitions, so a reader can skip the dense alignment.

## Open questions

- Dense mode: expand MAPLE input into a full alignment, or reject it and require sparse mode?
- Which MAPLE features are in scope: nucleotide data only, or also amino-acid data, which CMAPLE accepts?

## Related

- [kb/proposals/io-format-coverage.md](../proposals/io-format-coverage.md): ranking of format work
- [N-io-large-dataset-memory-constraint.md](N-io-large-dataset-memory-constraint.md): dense alignments increase peak memory
