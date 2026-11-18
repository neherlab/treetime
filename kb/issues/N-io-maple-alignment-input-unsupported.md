# MAPLE alignment input is not supported

TreeTime cannot read alignments in the MAPLE format. The format stores a reference sequence once and then the differences of each sample from that reference, so it is much smaller than FASTA for large sets of closely related genomes. CMAPLE reads it [[src](https://github.com/iqtree/cmaple/blob/3d45b1ab68e2d68a2825bf17a531e22200578cd6/alignment/alignment.h#L15-L18)], and IQ-TREE 3 can write it [[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/alignment/alignment.cpp#L6297)].

The per-sample difference lists map directly onto sparse partitions, so a reader can skip the dense alignment.

## Open questions

> [!IMPORTANT]
> **Decision required.** A MAPLE reader needs two scope choices that no decision covers:
>
> - Dense mode: expand MAPLE input into a full alignment, or reject MAPLE input in dense mode and require sparse mode
> - Alphabet: read nucleotide data only, or also amino-acid data, which CMAPLE accepts

## Related

- [kb/proposals/io-format-coverage.md](../proposals/io-format-coverage.md): ranking of format work
- [N-io-large-dataset-memory-constraint.md](N-io-large-dataset-memory-constraint.md): dense alignments increase peak memory

## Proposed implementation plan

The two scope choices above are open. Steps 1, 2, and 4 hold for every option; step 3 has one variant per dense-mode option.

1. Reader: add `packages/treetime-io/src/maple.rs` with a function that parses the `>REF` block and the per-sample difference lines (type, 1-based position, optional length) into a reference sequence plus a per-sample list of differences, with errors that name the line. Grammar oracle: CMAPLE `alignment/alignment.h` and its reader at the pinned revision above, cloned into `.repos/cmaple/`
2. Format selection: let `fn read_sequence_source()` (`packages/app-commands/src/commands/shared/sequence_inputs.rs`) pick the MAPLE reader by content (a first record named `REF` followed by difference lines) or by extension, beside `fasta_read_file`
3. Dense mode:
   - Expand: rebuild each sample's full sequence from the reference and its differences, then reuse the FASTA path unchanged
   - Reject: return an error that names sparse mode when a MAPLE file reaches a dense partition
4. Tests in `packages/treetime-io/src/__tests__/test_maple.rs`: a small MAPLE literal and the FASTA it encodes; the expanded sequences equal the FASTA records (oracle: the hand-written FASTA, cross-checked with IQ-TREE 3 conversion of the same file). Command-level check: sparse ancestral on a MAPLE file and on the equivalent FASTA give the same output
5. Validate: `just lint-rs`, `just test-rs`, and `just generated-check` when a CLI flag or schema changes
