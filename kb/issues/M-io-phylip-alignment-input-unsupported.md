# PHYLIP alignment input is not supported

The commands read alignments only as FASTA. A PHYLIP alignment fails in the FASTA reader:

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

## Decision needed

Implement a PHYLIP reader for parity with v0, or record the reduced input set as a decision in `kb/decisions/` and remove the bundled `aln.phylip`.
