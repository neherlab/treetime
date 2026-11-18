# `data/translate` maps every ambiguous codon to `X`

The script `data/translate`, which made the bundled peptide files (`data/*/translations/`), translates with Biopython's `forward_table.get(codon, "X")` (`data/translate`, line 105). The forward table holds only unambiguous codons, so every ambiguous codon becomes `X`, also when all its resolutions code one amino acid: `GCN` codes Ala, `ACN` codes Thr, `YTR` codes Leu. `fn _translate_codon()` (lines 100-105) maps `---` to `-`, the unambiguous stop codons to `*`, and every other codon through the forward table, so the ambiguous stop codon `TAR` also gives `X`.

The bundled peptides differ from the Nextclade translation rules at 464 tip residues of `rsv/a/2000` and 7,161 of `mpox/clade-ii/2000`. In every case the peptide file has `X` and the Nextclade rule gives an amino acid. Runs and tests that read these files see less tip information than a Nextclade translation of the same alignment gives.

## Fix direction

Translate with the Nextclade rules ([translate.rs#L8-L188](https://github.com/nextstrain/nextclade/blob/6190424cf8014d999fccf59d049c2c78a8bc28bb/packages/nextclade/src/translate/translate.rs#L8-L188)): an ambiguous codon whose resolutions all code one amino acid gives that amino acid; any other ambiguous or partly gapped codon gives `X`; `---` gives `-`. Regenerate the bundled peptide files, and update the tests and smoke snapshots that read them.

## Validation

`GCN`, `ACN`, `YTR` and `TAR` translate to `A`, `T`, `L` and `*`; `RAY` and `NNN` translate to `X`.

> [!IMPORTANT]
> **Investigation required.** The residue counts above (464 for `rsv/a/2000`, 7,161 for `mpox/clade-ii/2000`) were not reproduced. Recount them with the fixed script before and after regeneration.

## Proposed implementation plan

1. In `data/translate`, replace `fn _translate_codon()`: expand each IUPAC nucleotide code to its set of bases (Biopython `IUPACData.ambiguous_dna_values`), translate every resolution with the standard table (stop codons as `*`), and return the one amino acid when all resolutions agree; return `-` for `---` and `X` for every other codon with a gap, an `N`-only mix that disagrees, or disagreeing resolutions. This is the rule of Nextclade `translate.rs` cited above
2. Add a self-check to the script (or a small test beside it, following how `data/` scripts are tested, if any) with the cases of "Validation" (oracle: the standard genetic code and the Nextclade rule)
3. Regenerate `data/*/translations/` with the script, then recount the changed residues per dataset
4. Update the consumers of the peptide files: `dev/smoke.toml` snapshots (`just smoke` with a new baseline, reviewed case by case), `packages/treetime/src/test_utils/alignment.rs`, `packages/app-backend/src/__tests__/test_datasets.rs`, `packages/app-datasets/src/__tests__/test_examples.rs`; check `packages/app-cli/src/cli/command_help.rs` and the generated `packages/app-docs/src/generated/en/reference/cli.md` for any text that describes the translation
5. Validate: `just test-rs`, `just smoke`, and a diff of each regenerated file that shows only `X` changed to an amino acid or `*`
