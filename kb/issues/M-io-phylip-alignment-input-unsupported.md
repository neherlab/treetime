# PHYLIP alignment input is not supported

The commands read alignments only as FASTA, through `fn fasta_read_file()` [packages/treetime-io/src/fasta.rs#L13-L15](../../packages/treetime-io/src/fasta.rs#L13-L15). A PHYLIP alignment fails in the FASTA reader:

```
Error:
   0: FASTA input is incorrectly formatted: expected at least one FASTA record starting with character '>', but none found
```

v0 reads PHYLIP through `Bio.AlignIO`, and the bundled dataset `data/flu/h3n2/20` ships `aln.phylip`.

## Reproduction

```bash
treetime ancestral --tree=data/flu/h3n2/20/tree.nwk --alignment=data/flu/h3n2/20/aln.phylip --output-all=<dir>
```

This is the `ancestral/flu/h3n2/20/input-phylip-alignment` case of `dev/smoke`, declared an expected failure in `dev/smoke.toml`.

## Open question

> [!IMPORTANT]
> **Decision required.** v0 reads PHYLIP alignments and the bundled dataset ships one, but no decision covers the alignment formats of v1. Options:
>
> - Implement a PHYLIP reader, for parity with v0
> - Keep FASTA as the only alignment format: record the reduced input set in `kb/decisions/` and remove the bundled `aln.phylip` and its smoke case

## Proposed implementation plan

The format scope is open (Decision required above). One short plan per option. Note: a second PHYLIP file, `zika/86/zika.phylip.xz`, appears in the dataset test `packages/app-backend/src/__tests__/test_datasets.rs`.

- Option "implement a PHYLIP reader":
  1. Add `packages/treetime-io/src/phylip.rs` that reads the header (sample count, length) and the relaxed sequential and interleaved layouts into `Vec<FastaRecord>`, with errors that name the line and the count mismatch (format oracle: Biopython `Bio.AlignIO.PhylipIO`, the reader that v0 uses)
  2. Select the reader in `fn read_sequence_source()` (`packages/app-commands/src/commands/shared/sequence_inputs.rs`) by content (first line of two integers) or by the `.phylip`/`.phy` extension, beside `fasta_read_file()`
  3. Tests: sequential and interleaved literals of the same alignment give the records of the matching FASTA literal (oracle: Biopython `AlignIO.read(..., "phylip-relaxed")` on the same text)
  4. Remove the `fail =` line of the `input-phylip-alignment` case in `dev/smoke.toml`, and compare its output with the FASTA case of the same dataset (oracle: `aln.fasta.xz` of `data/flu/h3n2/20`)
- Option "FASTA only": write the decision file with user approval, delete `data/flu/h3n2/20/aln.phylip` and the smoke case, check the zika PHYLIP file in the dataset test, and make the FASTA error name PHYLIP when the first line is two integers
- Validate either option: `just lint-rs`, `just test-rs`, `just smoke`
