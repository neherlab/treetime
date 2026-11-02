# The augur node data of v1 writes indels on the amino-acid track and leaves them out on the nucleotide track

> [!IMPORTANT]
> **Decision required.** augur writes gap changes on both mutation tracks as positional tokens (`A123-`, `-123A`, `K100-`). v1 writes them on the amino-acid track and leaves them out on the nucleotide track. Choose one representation for both tracks; see "Decision required" below.

## augur behavior

augur writes every difference between the sequences of a parent and a child as one `<from><pos><to>` string, and treats the gap `-` as an ordinary aligned character on both tracks:

- **Nucleotides**: `def collect_mutations()` [[src](https://github.com/nextstrain/augur/blob/58f329042e7fc631ca275cf2df44c8720f80949f/augur/ancestral.py#L185)] writes the mutations of TreeTime v0 one string per site [[src](https://github.com/nextstrain/augur/blob/58f329042e7fc631ca275cf2df44c8720f80949f/augur/ancestral.py#L218-L219)]. v0 reconstructs with the `nuc` alphabet, in which the gap is a fifth state ([M-timetree-gap-state-v0-divergence.md](M-timetree-gap-state-v0-divergence.md)), so a deletion appears as `A123-` and a change from a gap as `-123A`
- **Amino acids**: a direct column comparison of the translations, without gap filtering, so a deleted residue appears as `K100-` ([kb/proposals/amino-acid-mutation-output-substitution-only.md](../proposals/amino-acid-mutation-output-substitution-only.md), section "Reference behavior")

augur export copies the node-data lists into the Auspice JSON, so deletions show on the branches in Auspice and can define clades.

## v1 behavior

- **Nucleotides**: `fn nucleotide_substitutions()` [packages/app-output/src/augur_node_data_ancestral.rs#L118-L137](../../packages/app-output/src/augur_node_data_ancestral.rs#L118-L137) keeps the substitutions of each edge and drops its insertion and deletion events. The written sequences keep their gaps
- **Amino acids**: `fn collect_aa_cds_node_data()` [packages/treetime/src/ancestral/aa.rs#L151](../../packages/treetime/src/ancestral/aa.rs#L151) chains substitutions with range-based indels, and `fn aa_mutation_strings()` [packages/app-output/src/augur_node_data_ancestral.rs#L139](../../packages/app-output/src/augur_node_data_ancestral.rs#L139) expands each range into one positional token per residue through `fn mutation_event_strings()` [packages/treetime/src/seq/mutation.rs#L216](../../packages/treetime/src/seq/mutation.rs#L216): `{state}{pos}-` for a deletion and `-{pos}{state}` for an insertion, the string form of augur. The substitution diff excludes gap positions, `fn is_reportable_sub()` [packages/treetime/src/ancestral/aa.rs#L224](../../packages/treetime/src/ancestral/aa.rs#L224), so gaps reach this track only through the indel events
- **Auspice JSON**: `fn group_mutations()` [packages/app-output/src/auspice.rs#L293](../../packages/app-output/src/auspice.rs#L293) drops the indels of the nucleotide track and keeps those of the amino-acid tracks

## Evidence

A four-tip alignment in which one tip has a three-base deletion:

```bash
printf '(((A:0.1,B:0.1):0.1,C:0.1):0.05,D:0.3);\n' > tree.nwk
printf '>A\nATGAAATTTCCCGGG\n>B\nATG---TTTCCCGGG\n>C\nATGAAATTTCCCGGA\n>D\nATGAAATTTCCAGGA\n' > aln.fasta
treetime ancestral --tree tree.nwk --alignment aln.fasta --method-anc marginal --gap-fill only-terminal --output-node-data v1.json
```

v1 writes `"muts": []` for `B`, where augur writes `A4-`, `A5-`, `A6-`. On the replayed Nextstrain production builds, gap changes are most of the nucleotide mutations of augur: 379 025 of the 386 075 mutations of the mpox B.1 build. The share of mutations that both engines report on the same node (Jaccard index) rises once gap changes are left out:

| Build | With gap changes | Without gap changes |
| --- | ---: | ---: |
| mpox clade I | 0.004 | 0.981 |
| mpox B.1 | 0.018 | 0.996 |
| seasonal-cov HKU1 | 0.297 | 0.806 |
| measles genome | 0.593 | 0.975 |
| ncov global, 6 months | 0.605 | 0.996 |
| hMPV A/G | 0.667 | 0.976 |

## Decision required

The options apply to both tracks together, so that one rule writes all mutation lists:

- **augur tokens on both tracks**: write the indel events of the nucleotide track through `fn mutation_event_strings()`, as the amino-acid track does. Matches augur for deletions; insertions (`-{pos}{state}`) need a check that their positions match augur's alignment columns
- **Substitutions only on both tracks**: drop the indel tokens of the amino-acid track as well, and record the divergence from augur in `kb/decisions/`. Deletions then disappear from Auspice and from clade definitions on both tracks
- **A separate indel output**: keep `muts` and `aa_muts` substitution-only and write indels in their own node-data field, which augur and Auspice do not read

The nucleotide model of the gap is a separate decision ([M-timetree-gap-state-v0-divergence.md](M-timetree-gap-state-v0-divergence.md)): with four states and indel ranges, v1 can still write augur tokens from its indel events. The coordinate question of amino-acid insertions is analyzed in [kb/proposals/amino-acid-mutation-output-substitution-only.md](../proposals/amino-acid-mutation-output-substitution-only.md).
