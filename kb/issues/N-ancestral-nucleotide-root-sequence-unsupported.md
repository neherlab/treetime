# Ancestral has no nucleotide root sequence

`augur ancestral --root-sequence <FASTA or GenBank>` reports the differences between a given nucleotide sequence, usually the reference genome, and the inferred root as mutations on the root branch, and translates it into the amino-acid root sequences when no `--aa-root-sequence` is given (`augur/ancestral.py`, `def _read_sequence_data()` and the `--root-sequence` help). TreeTime v1 `ancestral` has only the amino-acid counterpart, `--aa-root-sequence`: the root `aa_muts` of the augur node data list the differences between that reference and the inferred root. The root `muts` of the nucleotide sequence have no such reference.

## Who uses it

- nextstrain/seasonal-flu: `augur ancestral --root-sequence` with the FASTA reference of each segment (`workflow/snakemake_rules/core.smk`)
- nextstrain/rsv: `augur ancestral --root-sequence` with the GenBank reference of the build (`workflow/snakemake_rules/core.smk`)

A Nextstrain pipeline that runs `treetime ancestral` in place of `augur ancestral` cannot report the mutations of the root against the reference for these builds: the root of the node data has no nucleotide mutations.

## Fix direction

Add a nucleotide root sequence setting to `ancestral` next to `--aa-root-sequence`, with the same output rule for the root `muts` of the augur node data, and the translation of the nucleotide reference into the amino-acid root sequences when `--aa-root-sequence` is not given. A GenBank root sequence also needs a GenBank reader ([M-ancestral-genbank-annotation-unsupported.md](M-ancestral-genbank-annotation-unsupported.md)).

## Related issues

- [M-io-vcf-input-output-unimplemented.md](M-io-vcf-input-output-unimplemented.md): `--vcf-reference` plays the same role for VCF input in augur
