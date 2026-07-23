# Partition probability profiles use the opaque `dis` field name

`DenseSeqDistribution::dis` and `VarPos::dis` expose probability-profile arrays through an unexplained abbreviation [`packages/treetime/src/partition/dense.rs#L54`](../../packages/treetime/src/partition/dense.rs#L54) [`packages/treetime/src/partition/sparse.rs#L211`](../../packages/treetime/src/partition/sparse.rs#L211).

## Evidence

Call sites use forms such as `profile.dis`, `msg_to_parent.dis`, and `var.dis`. In the same modules, `profile` already names state-probability vectors and matrices. The abbreviation can also be read as distance, dissimilarity, or distribution, so its meaning depends on the surrounding type.

The field is part of central dense and sparse partition structures. A local rename would leave two terms for the same scientific object and make conversions harder to audit.

## Required vocabulary

Use `profile` for state-probability arrays throughout dense and sparse partition APIs. Keep `distribution` for time or branch-length distribution objects, and keep distance names explicit where they represent a metric.

## Validation

- Semantic reference search confirms every `.dis` access is classified before renaming.
- Dense and sparse reconstruction tests compare complete partition values before and after the rename.
- Public serialized field names, if any, are checked explicitly; v1 has not shipped, so no compatibility alias is needed.
