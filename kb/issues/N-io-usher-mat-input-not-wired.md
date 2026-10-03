# UShER MAT input is not wired to any command

`fn usher_mat_pb_read()` [packages/util-usher-mat/src/lib.rs#L56](../../packages/util-usher-mat/src/lib.rs#L56) parses MAT protobuf, but no command reads a MAT. The converter from MAT to the graph was removed with other unreachable read paths in commit `b944399c`. A user who has a MAT, for example one of the daily UCSC trees, must convert it to Newick and FASTA before TreeTime can date it, and that conversion expands the sparse mutation data into a dense alignment.

## Required behavior for a future reader

- Keep the distinction between an absent branch length and an explicit zero in the embedded Newick tree. The removed converter turned an absent length into `0.0`, so unknown distance looked like observed zero change
- Build sparse partitions from the per-node mutations and a reference sequence, without a dense alignment

## Open questions

- Where does the reference sequence come from? A MAT stores reference states only at mutated positions, so the full sequence needs a separate input
- How are masked sites and missing data in the MAT represented in the partition?
- Which commands accept MAT input: `timetree` and `clock` only, or every tree-reading command?

## Related

- [kb/proposals/io-format-coverage.md](../proposals/io-format-coverage.md): ranking of format work
- [kb/proposals/unified-input-format-support.md](../proposals/unified-input-format-support.md): input paths that build partitions from formats with embedded data
- [kb/decisions/io-usher-mat-gaps-as-missing-data.md](../decisions/io-usher-mat-gaps-as-missing-data.md): how the MAT writers handle gaps
