# Amino-acid mutation output indel representation is undecided

> [!IMPORTANT]
> **Decision required.** Restricting the amino-acid track to substitutions drops the deletion tokens that augur also writes as positional gap-substitutions, and keeping range-based indels requires verifying their coordinates against augur's alignment columns. See "Open question" and the linked proposal.

The nucleotide mutation output is substitution-only, matching augur node-data. The amino-acid track still emits range-based insertion and deletion events, expanded to per-residue tokens. Whether the amino-acid track should be restricted to substitutions (as the nucleotide track is) is an open design question, not yet decided by the team.

## Evidence

- Amino-acid node-data chains substitutions with range-based indels: `fn collect_aa_cds_node_data()` [`packages/treetime/src/ancestral/aa.rs#L158-L204`](../../packages/treetime/src/ancestral/aa.rs#L158-L204) and `fn aa_mutation_strings()` [`packages/app-output/src/augur_node_data_ancestral.rs#L139-L147`](../../packages/app-output/src/augur_node_data_ancestral.rs#L139-L147). The substitution diff excludes gap positions, `fn is_reportable_sub()` [`packages/treetime/src/ancestral/aa.rs#L231-L234`](../../packages/treetime/src/ancestral/aa.rs#L231-L234), so gaps reach output only via the indel track.
- The Auspice writer drops indels for the nucleotide track but not amino-acid tracks: `fn group_mutations()` [`packages/app-output/src/auspice.rs#L266-L280`](../../packages/app-output/src/auspice.rs#L266-L280).
- `fn mutation_event_strings()` [`packages/treetime/src/seq/mutation.rs#L217-L233`](../../packages/treetime/src/seq/mutation.rs#L217-L233) expands a deletion into one `{state}{pos}-` token per residue (1-based), which is augur's positional gap-substitution form (`K100-`). An insertion becomes `-{pos}{state}` tokens.

## Open question

augur represents amino-acid deletions as positional gap-substitutions (`K100-`), not as range events. TreeTime deletion tokens already use the same string form, so restricting the amino-acid track to substitutions would drop tokens that augur writes and diverge from augur. Keeping range-based indels requires verifying that their coordinates match augur's alignment columns, especially for insertions, which have no augur counterpart. The decision and its options are analyzed in [kb/proposals/amino-acid-mutation-output-substitution-only.md](../proposals/amino-acid-mutation-output-substitution-only.md).

Implementation waits until the team selects an option.
