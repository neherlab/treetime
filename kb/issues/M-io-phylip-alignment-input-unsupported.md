# PHYLIP alignment input is not supported

The commands read alignments only as FASTA, through `fn fasta_read_file()` [packages/treetime-io/src/fasta.rs#L13-L15](../../packages/treetime-io/src/fasta.rs#L13-L15). A PHYLIP alignment fails in the FASTA reader:

```
Error:
   0: FASTA input is incorrectly formatted: expected at least one FASTA record starting with character '>', but none found
```

v0 reads PHYLIP through `Bio.AlignIO`, and the bundled dataset `data/flu/h3n2/20` ships `aln.phylip`.

## Reproduction

```bash
treetime ancestral --tree=data/flu/h3n2/20/tree.nwk --aln=data/flu/h3n2/20/aln.phylip --output-all=<dir>
```

This is the `ancestral/flu/h3n2/20/input-phylip-alignment` case of `dev/smoke`, declared an expected failure in `dev/smoke.toml`.

## Open question

> [!IMPORTANT]
> **Decision required.** v0 reads PHYLIP alignments and the bundled dataset ships one, but no decision covers the alignment formats of v1. Options:
>
> - Implement a PHYLIP reader, for parity with v0
> - Keep FASTA as the only alignment format: record the reduced input set in `kb/decisions/` and remove the bundled `aln.phylip` and its smoke case
