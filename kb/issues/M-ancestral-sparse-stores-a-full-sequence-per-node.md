# Sparse reconstruction stores a full sequence per node

## Symptom

The sparse backend materializes one full-length `Seq` for every node in `SparseSeqInfo::sequence`, plus a second one for every tip in `SparseNodePartition::emitted`. Storage is `O(N * L)` in node count and alignment length, which is the cost the sparse representation exists to avoid for probability vectors.

On `data/sc2/4500` (10199 nodes, L = 29903): 305 MB of per-node sequences and 153 MB of tip `emitted`, against 0.21 MB for all 26300 edge substitutions and 1.95 MB for all missing-data ranges. Roughly 458 MB where 2.2 MB carries the same information.

`capture_ancestral_states` ([timetree/convergence/sequence_changes.rs#L45](../../packages/treetime/src/timetree/convergence/sequence_changes.rs#L45)) additionally snapshots every internal sequence twice per timetree iteration, about 610 MB on the same dataset, purely to count changed positions.

## Reproduction

Run `ancestral` on `data/sc2/4500` and measure peak RSS, or compute the accounting directly: node count times alignment length against the sum of the per-edge mutation lists and per-node range tracks.

## Impact and scope

Memory, not correctness. It bounds the alignment length and tree size the sparse path can handle, and the ratio worsens with alignment length at fixed divergence - exactly the regime the sparse backend targets. Dense is unaffected, since its `(L, K)` profiles dominate regardless.

## Root cause

The stored sequences are redundant with data already held elsewhere. A node's sequence is determined by the root sequence, the substitutions and indels on the path to it, and its own missing-data mask.

This was verified per edge on `data/sc2/4500`: across all 10198 edges there are no reported substitutions that are not actual differences, and no positions differing between parent and child with both endpoints canonical that lack a reported substitution. Every unreported difference has `N`, `-`, or an IUPAC code at one end, and each of those is carried by a track that is already stored. The "parent masked, child canonical" case that would lose a substitution is unreachable, because an internal node's `non_char` is the intersection of its children's.

## Fix approach

See the proposal [Mutation-first sequence representation](../proposals/mutation-first-sequence-representation.md) for the design axes, staging and validation plan. The open axes are the random-access `PartitionBranchOps::node_sequence` accessor, `--sample-from-profile=all` (whose realization is not compressible), and whether the augur per-node `sequence` field can be gated behind a flag.

The first stage is independent and has no open questions: make the timetree convergence check compare per-node mutation sets rather than sequences.
