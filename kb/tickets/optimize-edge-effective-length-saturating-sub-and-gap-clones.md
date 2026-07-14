# edge_effective_length saturating_sub clamps to zero and clones gap vectors

Two related issues in edge_effective_length.

## Instances

### saturating_sub clamps to 0

v1: [`packages/treetime/src/partition/marginal_sparse.rs#L255`](../../packages/treetime/src/partition/marginal_sparse.rs#L255) and [`marginal_dense.rs#L153`](../../packages/treetime/src/partition/marginal_dense.rs#L153)

Returns 0 when non-char positions exceed sequence length. Used as denominator in `edge_subs().len() / edge_effective_length()`, producing division by zero.

### range_union clones gap vectors

v1: [`packages/treetime/src/partition/marginal_dense.rs#L122`](../../packages/treetime/src/partition/marginal_dense.rs#L122) and [`marginal_sparse.rs#L380`](../../packages/treetime/src/partition/marginal_sparse.rs#L380)

`fn range_union` clones both gap vectors on every call. Allocation pressure on heavily gapped alignments.

## Related issues

- Source: [kb/issues/N-code-quality-conventions.md](../issues/N-code-quality-conventions.md) -- delete after full resolution
