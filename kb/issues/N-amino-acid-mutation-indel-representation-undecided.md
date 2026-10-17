# Amino-acid mutation output indel representation is undecided

> [!IMPORTANT]
> **Decision required.** Restricting the amino-acid track to substitutions drops augur's positional gap-substitution tokens, and keeping range-based indels requires verifying their coordinates against augur's alignment columns. See "Open question" and the linked proposal.

The nucleotide mutation output is substitution-only, matching augur node-data. The amino-acid track still emits range-based insertion and deletion events, expanded to per-residue tokens. Whether the amino-acid track should be restricted to substitutions (as the nucleotide track is) is an open design question, not yet decided by the team.

## Evidence

- Amino-acid node-data chains substitutions with range-based indels [`packages/treetime/src/ancestral/aa.rs#L153-L199`](../../packages/treetime/src/ancestral/aa.rs#L153-L199) and [`packages/app-output/src/augur_node_data_ancestral.rs#L139-L147`](../../packages/app-output/src/augur_node_data_ancestral.rs#L139-L147); the substitution diff excludes gap positions [`packages/treetime/src/ancestral/aa.rs#L226-L229`](../../packages/treetime/src/ancestral/aa.rs#L226-L229), so gaps reach output only via the indel track.
- The Auspice writer drops indels for the nucleotide track but not amino-acid tracks `fn group_mutations()` [`packages/app-output/src/auspice.rs#L266-L280`](../../packages/app-output/src/auspice.rs#L266-L280).

## Open question

augur represents amino-acid deletions as positional gap-substitutions (`K100-`), not as range events, so restricting the amino-acid track to substitutions would drop those tokens and diverge from augur; keeping range-based indels requires verifying that their coordinates match augur's alignment columns, especially for insertions. The decision and its options are analyzed in [kb/proposals/amino-acid-mutation-output-substitution-only.md](../proposals/amino-acid-mutation-output-substitution-only.md).

Implementation waits until the team selects an option.
