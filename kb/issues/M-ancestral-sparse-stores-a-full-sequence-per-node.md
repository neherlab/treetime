# Sparse reconstruction stores a full sequence per node

## Symptom

The sparse backend materializes one full-length `Seq` for every node in `SparseNodeState::sequence` ([packages/treetime/src/partition/storage/sparse.rs](../../packages/treetime/src/partition/storage/sparse.rs)), and during partition construction in `FitchSeqInfo::sequence`. Storage is `O(N * L)` in node count and alignment length, which is the cost the sparse representation exists to avoid for probability vectors.

On `data/sc2/4500` (10199 nodes, L = 29903): 305 MB of per-node sequences, against 0.21 MB for all 26300 edge substitutions and 1.95 MB for all missing-data ranges.

The marginal forward pass keeps a node's sequence when it rebuilds it unchanged, and the timetree convergence snapshot (`fn capture_ancestral_states()` in [packages/treetime/src/timetree/convergence/sequence_changes.rs](../../packages/treetime/src/timetree/convergence/sequence_changes.rs)) shares the node sequences and stores only the states inferred at variable positions, so neither adds a second set of sequences.

## Reproduction

Run `ancestral` on `data/sc2/4500` and measure peak RSS, or compute the accounting directly: node count times alignment length against the sum of the per-edge mutation lists and per-node range tracks.

## Impact and scope

Memory, not correctness. It bounds the alignment length and tree size the sparse path can handle, and the ratio worsens with alignment length at fixed divergence - exactly the regime the sparse backend targets. Dense is unaffected, since its `(L, K)` profiles dominate regardless.

## Root cause

The stored sequences are redundant with data already held elsewhere. A node's sequence is determined by the root sequence, the substitutions and indels on the path to it, and its own missing-data mask.

This was verified per edge on `data/sc2/4500`: across all 10198 edges there are no reported substitutions that are not actual differences, and no positions differing between parent and child with both endpoints canonical that lack a reported substitution. Every unreported difference has `N`, `-`, or an IUPAC code at one end, and each of those is carried by a track that is already stored. The "parent masked, child canonical" case that would lose a substitution is unreachable, because an internal node's `non_char` is the intersection of its children's.

## Fix approach

See the proposal [Mutation-first sequence representation](../proposals/mutation-first-sequence-representation.md) for the design axes, staging and validation plan. The open axes are the random-access `PartitionBranchOps::node_sequence` accessor, `--sample-from-profile=all` (whose realization is not compressible), and whether the augur per-node `sequence` field can be gated behind a flag.
