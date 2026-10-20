# Partition sequence state is exposed as `seq.sequence`

`DenseNodeState.seq` and `FitchNodeData.seq` hold `DenseSeqInfo` and `FitchSeqInfo`, whose actual character storage is named `sequence` [packages/treetime/src/partition/storage/dense.rs#L10-L14](../../packages/treetime/src/partition/storage/dense.rs#L10-L14), [packages/treetime/src/partition/storage/dense.rs#L58-L65](../../packages/treetime/src/partition/storage/dense.rs#L58-L65), [packages/treetime/src/partition/storage/sparse.rs#L225-L228](../../packages/treetime/src/partition/storage/sparse.rs#L225-L228), [packages/treetime/src/partition/storage/sparse.rs#L269-L277](../../packages/treetime/src/partition/storage/sparse.rs#L269-L277). The sparse marginal state `SparseNodeState` holds `sequence` directly and does not have the problem [packages/treetime/src/partition/storage/sparse.rs#L65-L69](../../packages/treetime/src/partition/storage/sparse.rs#L65-L69).

## Evidence

Call sites consequently use `node.seq.sequence`: 14 uses in 5 files, in the Fitch partition and passes, the dense marginal partition, the shared marginal pass, and the Fitch tests. The neighboring indel-interval fields (`gaps`, `unknown`, `non_char`, `variable_indel`) and, in `FitchSeqInfo`, the composition and the Fitch state distribution show that the outer field represents a bundle of sequence reconstruction state rather than a `Seq` value. The short name creates both stuttering and a false type cue.

Readers repeatedly have to infer which layer of sequence state is being accessed.

## Required vocabulary

Rename the outer field for its complete role, such as `sequence_state`, while retaining `sequence` for the actual character sequence. Dense and Fitch structures must use the same concept name.

## Validation

- Semantic impact analysis covers field construction, destructuring, and macro-generated access.
- Dense and Fitch whole-structure tests confirm a naming-only change.
- Documentation distinguishes alignment input, reconstructed sequence, probability profile, and per-node sequence state.
