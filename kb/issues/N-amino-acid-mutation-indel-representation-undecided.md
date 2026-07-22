# Amino-acid mutation output indel representation is undecided

The nucleotide mutation output is substitution-only, matching augur node-data. The amino-acid track still emits range-based insertion and deletion events, expanded to per-residue tokens. Whether the amino-acid track should be restricted to substitutions (as the nucleotide track is) is an open design question, not yet decided by the team.

## Evidence

- Amino-acid node-data chains substitutions with range-based indels [`packages/treetime/src/commands/ancestral/aa_node_data.rs#L270-L276`](../../packages/treetime/src/commands/ancestral/aa_node_data.rs#L270-L276); the substitution diff excludes gap positions [`packages/treetime/src/commands/ancestral/aa_node_data.rs#L326-L329`](../../packages/treetime/src/commands/ancestral/aa_node_data.rs#L326-L329), so gaps reach output only via the indel track.
- The Auspice writer drops indels for the nucleotide track but not amino-acid tracks `fn group_mutations()` [`packages/treetime/src/commands/shared/tree_output.rs#L1473-L1502`](../../packages/treetime/src/commands/shared/tree_output.rs#L1473-L1502).

## Open question

augur represents amino-acid deletions as positional gap-substitutions (`K100-`), not as range events, so restricting the amino-acid track to substitutions would drop those tokens and diverge from augur; keeping range-based indels requires verifying that their coordinates match augur's alignment columns, especially for insertions. The decision and its options are analyzed in [kb/proposals/amino-acid-mutation-output-substitution-only.md](../proposals/amino-acid-mutation-output-substitution-only.md).

Blocked on a team decision; no ticket until an option is selected.
