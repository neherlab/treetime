# `--aa-model jtt92` writes residues at the stop codon of internal nodes

With `--aa-model jtt92`, the amino-acid route writes an amino acid at the annotated stop codon of every internal node, while every tip has `X` there. The residue is `M` in F, OPG001 and OPG105, and `Q` in L. This affects 1,997 to 1,999 internal nodes per CDS on `rsv/a/2000` and 1,987 per CDS on `mpox/clade-ii/2000`. With `jtt92` the route also reconstructs none of the internal stop codons (57 on `rsv/a/2000`, 142 on `mpox/clade-ii/2000`) that the default model reconstructs.

`jtt92` uses the 20-state alphabet without a stop state (`AlphabetName::AaNoStop`, [packages/app-commands/src/commands/ancestral/aa_model.rs#L32-L34](../../packages/app-commands/src/commands/ancestral/aa_model.rs#L32-L34)), so every stop codon of the tips becomes `X`.

## Reproduction

```
treetime ancestral --tree data/rsv/a/2000/tree.nwk --alignment data/rsv/a/2000/aln.fasta.xz \
  --cds-selection NS1,NS2,N,P,M,SH,G,F,M2-1,M2-2,L \
  --translations 'data/rsv/a/2000/translations/{cds}.fasta.xz' --aa-model jtt92 \
  --output-translations '<dir>/{cds}.fasta'
```

The last residue of each peptide is `X` at every tip and `M` or `Q` at every internal node.

> [!IMPORTANT]
> **Investigation required.** The mechanism is not traced. A probable cause: the terminal stop column has `X` at every tip, so the sparse reconstruction treats it as a column without data and writes the most probable state of the equilibrium frequencies at internal nodes (`fn GloballyMissing.new()` [packages/treetime/src/partition/marginal/sparse/globally_missing.rs#L17-L23](../../packages/treetime/src/partition/marginal/sparse/globally_missing.rs#L17-L23)). This does not explain why the residue differs between CDS. Confirm the mechanism, and check whether internal stop codons of other columns become a residue in the same way.

> [!IMPORTANT]
> **Decision required.** If the cause is the fill of columns without data, the same rule writes the most probable equilibrium state at internal nodes for every column that has no data at any tip, in every alphabet. Options: write the unknown character (`X`, `N`) at internal nodes for such columns in every alphabet; keep the fill and exclude only stop-codon columns from 20-state models; or reject `jtt92` for CDS whose tips carry stop codons.

## Related

- [M-ancestral-aa-outputs-disagree-with-translated-nucleotide-outputs.md](M-ancestral-aa-outputs-disagree-with-translated-nucleotide-outputs.md): the recommended design translates TreeTime's own codons, so empirical amino-acid models remain only for protein-only input
