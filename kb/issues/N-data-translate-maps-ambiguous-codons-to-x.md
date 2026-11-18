# `data/translate` maps every ambiguous codon to `X`

The script `data/translate`, which made the bundled peptide files (`data/*/translations/`), translates with Biopython's `forward_table.get(codon, "X")` (`data/translate`, line 105). The forward table holds only unambiguous codons, so every ambiguous codon becomes `X`, also when all its resolutions code one amino acid: `GCN` codes Ala, `ACN` codes Thr, `YTR` codes Leu.

The bundled peptides differ from the Nextclade translation rules at 464 tip residues of `rsv/a/2000` and 7,161 of `mpox/clade-ii/2000`. In every case the peptide file has `X` and the Nextclade rule gives an amino acid. Runs and tests that read these files see less tip information than a Nextclade translation of the same alignment gives.

## Fix direction

Translate with the Nextclade rules ([translate.rs#L8-L188](https://github.com/nextstrain/nextclade/blob/6190424cf8014d999fccf59d049c2c78a8bc28bb/packages/nextclade/src/translate/translate.rs#L8-L188)): an ambiguous codon whose resolutions all code one amino acid gives that amino acid; any other ambiguous or partly gapped codon gives `X`; `---` gives `-`. Regenerate the bundled peptide files, and update the tests and smoke snapshots that read them.

## Validation

`GCN`, `ACN`, `YTR` and `TAR` translate to `A`, `T`, `L` and `*`; `RAY` and `NNN` translate to `X`.
