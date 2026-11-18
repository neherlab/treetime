# UShER MAT input is not wired to any command

`fn usher_mat_pb_read()` [packages/util-usher-mat/src/lib.rs#L56](../../packages/util-usher-mat/src/lib.rs#L56) parses MAT protobuf, but no command reads a MAT. The converter from MAT to the graph was removed with other unreachable read paths in commit `b944399c`. A user who has a MAT, for example one of the daily UCSC trees, must convert it to Newick and FASTA before TreeTime can date it, and that conversion expands the sparse mutation data into a dense alignment.

## Required behavior for a future reader

- Keep the distinction between an absent branch length and an explicit zero in the embedded Newick tree. The removed converter turned an absent length into `0.0`, so unknown distance looked like observed zero change
- Build sparse partitions from the per-node mutations and a reference sequence, without a dense alignment

## Open questions

> [!IMPORTANT]
> **Decision required.** The input contract of a MAT reader is not approved. Open points:
>
> - Source of the reference sequence. A MAT stores reference states only at mutated positions, so the full sequence needs a separate input (for example a FASTA file of the reference)
> - Representation of masked sites and missing data of the MAT in the partition
> - Commands that accept MAT input: `timetree` and `clock` only, or every command that reads a tree

## Related

- [kb/proposals/io-format-coverage.md](../proposals/io-format-coverage.md): ranking of format work
- [kb/proposals/unified-input-format-support.md](../proposals/unified-input-format-support.md): input paths that build partitions from formats with embedded data
- [kb/decisions/io-usher-mat-gaps-as-missing-data.md](../decisions/io-usher-mat-gaps-as-missing-data.md): how the MAT writers handle gaps

## Proposed implementation plan

The steps below stop at the open decision above; steps 2 to 4 need its answers.

1. Restore a MAT-to-graph reader in `packages/treetime-io/src/usher_mat.rs` (the writer lives there; commit `b944399c` removed the old reader from this file and from `packages/treetime/src/commands/timetree/initialization.rs`). Parse the embedded Newick with the project Newick reader so that an absent branch length stays `None`, and attach the per-node mutations of `UsherTree` (`usher_mat_pb_read()` in `packages/util-usher-mat/src/lib.rs`) to the graph edges
2. Build sparse partitions from the per-node mutations and the reference sequence input chosen in the decision, with the chosen representation of masked sites and missing data (see `kb/decisions/io-usher-mat-gaps-as-missing-data.md` for the writer side)
3. Add a MAT input option to the commands chosen in the decision, through the shared input path described in `kb/proposals/unified-input-format-support.md`; regenerate schemas with `just gen` and check with `just generated-check`
4. Tests in `packages/treetime-io/src/__tests__/test_usher_mat.rs`:
   - Round trip: write a small graph with sparse partitions through the MAT writer, read it back, compare the whole graph and the mutation sets (oracle: the input graph)
   - Branch lengths: a Newick string with one absent length and one explicit `0.0` keeps `None` and `Some(0.0)` (oracle: the Newick text)
   - Sparse equivalence: a small hand-written MAT plus reference gives the same leaf sequences as the equivalent FASTA alignment through the dense path (oracle: dense vs sparse comparison)
5. Validate with `just lint-rs`, `just test-rs`, and `just r treetime timetree` on a MAT made from a `data/` dataset, compared with the Newick and FASTA run of the same data
