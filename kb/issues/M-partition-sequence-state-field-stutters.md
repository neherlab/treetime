# Partition sequence state is exposed as `seq.sequence`

`DenseNodePartition.seq` and `SparseNodePartition.seq` hold `DenseSeqInfo` and `SparseSeqInfo`, whose actual character storage is named `sequence` [`packages/treetime/src/partition/dense.rs#L10`](../../packages/treetime/src/partition/dense.rs#L10) [`packages/treetime/src/partition/sparse.rs#L16`](../../packages/treetime/src/partition/sparse.rs#L16).

## Evidence

Call sites consequently use `node.seq.sequence`. The neighboring indel-interval field and profile-related state show that the outer field represents a bundle of sequence reconstruction state rather than a `Seq` value. The short name creates both stuttering and a false type cue.

This naming appears in marginal passes, reconstruction, output projection, and tests, so readers repeatedly have to infer which layer of sequence state is being accessed.

## Required vocabulary

Rename the outer field for its complete role, such as `sequence_state`, while retaining `sequence` for the actual character sequence. Dense and sparse structures must use the same concept name.

## Validation

- Semantic impact analysis covers field construction, destructuring, and macro-generated access.
- Dense and sparse whole-structure tests confirm a naming-only change.
- Documentation distinguishes alignment input, reconstructed sequence, probability profile, and per-node sequence state.
