# UShER MAT output fails on every run with amino-acid reconstruction

`ancestral` with `--translations` passes the amino-acid mutations of every node to the MAT writers (`mat-pb`, `mat-json`) together with the nucleotide mutations [packages/app-output/src/ancestral_tree_output.rs#L198](../../packages/app-output/src/ancestral_tree_output.rs#L198). `fn mat_mutation()` rejects every mutation that is not on the nucleotide track [packages/app-output/src/tree_output.rs#L515-L517](../../packages/app-output/src/tree_output.rs#L515-L517), so the whole command stops:

```
treetime ancestral --tree=data/rsv/a/20/tree.nwk --aln=data/rsv/a/20/aln.fasta.xz --method-anc=marginal \
  --translations='data/rsv/a/20/translations/%GENE.fasta.xz' --cdses=NS1,NS2 \
  --output-all=<dir> --output-selection=mat-json

Error:
   0: Node 'PP377517' has an amino-acid mutation that UShER MAT cannot represent
```

With `--output-selection=all`, the command exits with this error after it has written the outputs that come before MAT, so the MAT files, `graph.json` and the Graphviz file are missing.

## Reference behavior

A MAT stores nucleotide mutations only: the `mut` message of `parsimony.proto` has a position, a chromosome and nucleotide codes (0 to 3 for A, C, G, T), and no protein track ([UShER](https://github.com/yatisht/usher), revision `ac9c982d`). matUtils derives amino-acid changes from the nucleotide mutations and a GTF annotation (`matUtils summary --translate`).

## Open question

How should the MAT writers handle amino-acid mutations?

- Keep the error, and document that MAT output cannot be combined with `--translations`
- Write the nucleotide mutations only and leave the amino-acid track out, with or without a warning. The amino-acid changes stay in Auspice JSON, augur node-data JSON and the reconstructed amino-acid FASTA

## Smoke coverage

The `translations-mat` row of `dev/smoke.toml` runs MAT output with `--translations` and declares this failure. The other amino-acid rows leave MAT out of their selection.

## Related

- [kb/decisions/io-usher-mat-gaps-as-missing-data.md](../decisions/io-usher-mat-gaps-as-missing-data.md): how the MAT writers handle insertions and deletions
