# Partition probability profiles use the opaque `dis` field name

Three partition structures expose probability-profile arrays through the unexplained abbreviation `dis`:

- `DenseSeqDistribution::dis` [`packages/treetime/src/partition/storage/dense.rs#L85`](../../packages/treetime/src/partition/storage/dense.rs#L85)
- `FixedProfile::dis` [`packages/treetime/src/partition/storage/sparse.rs#L220`](../../packages/treetime/src/partition/storage/sparse.rs#L220), with the method `FixedProfiles::map_dis()` [`packages/treetime/src/partition/storage/sparse.rs#L189`](../../packages/treetime/src/partition/storage/sparse.rs#L189)
- `VarPos::dis` [`packages/treetime/src/partition/storage/sparse.rs#L264`](../../packages/treetime/src/partition/storage/sparse.rs#L264)

## Evidence

Accesses to `.dis` occur in 24 files of `packages/treetime/src`. Call sites use forms such as `profile.dis`, `var_pos.dis`, `fixed.dis`, and `msg_to_parent.dis`. In the same modules, `profile` already names state-probability vectors and matrices. The abbreviation can also be read as distance, dissimilarity, or distribution, so its meaning depends on the surrounding type.

The field is part of central dense and sparse partition structures. A local rename would leave two terms for the same scientific object and make conversions harder to audit.

## Required vocabulary

Use `profile` for state-probability arrays throughout dense and sparse partition APIs. Keep `distribution` for time or branch-length distribution objects, and keep distance names explicit where they represent a metric.

## Validation

- Semantic reference search confirms every `.dis` access and every `dis` parameter or method name (`map_dis()`, `VarPosMap::insert()`, `FixedProfiles::insert()`) is classified before renaming.
- Dense and sparse reconstruction tests compare complete partition values before and after the rename.
- Public serialized field names, if any, are checked explicitly; v1 has not shipped, so no compatibility alias is needed.
