# UShER MAT output fails on every run with amino-acid reconstruction

`ancestral` with `--translations` passes the amino-acid mutations of every node to the MAT writers (`mat-pb`, `mat-json`) together with the nucleotide mutations [packages/app-output/src/ancestral_tree_output.rs#L198](../../packages/app-output/src/ancestral_tree_output.rs#L198). `fn mat_mutation()` rejects every mutation that is not on the nucleotide track [packages/app-output/src/tree_output.rs#L511-L513](../../packages/app-output/src/tree_output.rs#L511-L513), so the whole command stops:

```
treetime ancestral --tree=data/rsv/a/20/tree.nwk --aln=data/rsv/a/20/aln.fasta.xz --method-anc=marginal \
  --translations='data/rsv/a/20/translations/%GENE.fasta.xz' --cdses=NS1,NS2 \
  --output-all=<dir> --output-selection=mat-json

Error:
   0: Node 'PP377517' has an amino-acid mutation that UShER MAT cannot represent
```

With `--output-selection=all`, the command exits with this error after it has written the outputs that come before MAT, so the MAT files, `graph.json` and the Graphviz file are missing.

## Reference behavior

A MAT stores nucleotide mutations only: the `mut` message of `parsimony.proto` has a position, a chromosome and nucleotide codes (0 to 3 for A, C, G, T), and no protein track [[src](https://github.com/yatisht/usher/blob/ac9c982d937c3bc43e2c16a0b73d30cf2b937118/parsimony.proto#L4-L11)]. The schema marks the `chromosome` string as unused. The UShER reader turns each state into a 4-bit nucleotide mask and asserts that every state is between 0 and 3 [[src](https://github.com/yatisht/usher/blob/ac9c982d937c3bc43e2c16a0b73d30cf2b937118/src/mutation_annotated_tree.cpp#L77-L85)], so no encoding of the 20 amino-acid states fits the format.

UShER tools compute amino-acid changes when they read a MAT: `matUtils summary --translate` and `matUtils extract --write-taxodium` both take the MAT, a GTF annotation and the reference FASTA [[src](https://github.com/yatisht/usher/blob/ac9c982d937c3bc43e2c16a0b73d30cf2b937118/src/matUtils/summary.cpp#L17-L31)] [[src](https://github.com/yatisht/usher/blob/ac9c982d937c3bc43e2c16a0b73d30cf2b937118/src/matUtils/extract.cpp#L13-L14)]. Writing the nucleotide mutations therefore loses no amino-acid information that a MAT could hold.

Taxonium JSONL stores amino-acid and nucleotide mutations together, and is the only format of the UShER family that can carry TreeTime's amino-acid track: [N-io-taxonium-jsonl-output-unsupported.md](N-io-taxonium-jsonl-output-unsupported.md).

## Open question

How should the MAT writers handle amino-acid mutations?

- Keep the error, and document that MAT output cannot be combined with `--translations`
- Write the nucleotide mutations only and leave the amino-acid track out, with or without a warning. The amino-acid changes stay in Auspice JSON, augur node-data JSON and the reconstructed amino-acid FASTA

Recommendation, not yet approved: write the nucleotide mutations only, with a warning. The error stops MAT output for every run with `--translations`, although the MAT could never hold the amino-acid track, and UShER tools compute the amino-acid changes again from the nucleotide mutations.

## Smoke coverage

The `translations-mat` row of `dev/smoke.toml` runs MAT output with `--translations` and declares this failure. The other amino-acid rows leave MAT out of their selection.

## Related

- [kb/decisions/io-usher-mat-gaps-as-missing-data.md](../decisions/io-usher-mat-gaps-as-missing-data.md): how the MAT writers handle insertions and deletions
