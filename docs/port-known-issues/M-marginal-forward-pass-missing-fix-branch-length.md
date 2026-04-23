# Marginal forward pass missing fix_branch_length

The backward marginal pass applies `fix_branch_length()` to clamp zero-length branches to a minimum value, but the forward pass reads the raw branch length from the graph payload without clamping. This means backward and forward passes evaluate the same edge under different transition matrices when the raw branch length is below the minimum threshold.

## Affected code

Sparse marginal:

- Backward: [packages/treetime/src/representation/partition/marginal_passes.rs#L135-L136](../../packages/treetime/src/representation/partition/marginal_passes.rs#L135-L136) applies `fix_branch_length()`
- Forward: [packages/treetime/src/representation/partition/marginal_passes.rs#L202](../../packages/treetime/src/representation/partition/marginal_passes.rs#L202) uses raw `branch_length`

Dense marginal:

- Backward: [packages/treetime/src/representation/partition/marginal_dense.rs#L278-L279](../../packages/treetime/src/representation/partition/marginal_dense.rs#L278-L279) applies `fix_branch_length()`
- Forward: [packages/treetime/src/representation/partition/marginal_dense.rs#L300](../../packages/treetime/src/representation/partition/marginal_dense.rs#L300) uses raw `branch_length`

## Impact

For edges with branch length exactly zero, the backward pass uses `fix_branch_length(length, 0.0)` which produces a small positive value, while the forward pass uses `0.0`. The forward-pass transition matrix `exp(Q * 0) = I` (identity), meaning the parent message passes through unchanged. The backward-pass transition matrix `exp(Q * epsilon)` is close to identity but not exact. The resulting posterior combines messages computed under two different models.

The practical impact is small because `fix_branch_length` produces a very small minimum value (a fraction of one mutation per sequence length). The transition matrices differ by ~1e-6 or less.

## Fix

Apply `fix_branch_length()` consistently in both passes, or remove it from the backward pass. The choice depends on whether zero-length branches should be treated as literally zero (identity matrix) or as a small positive value (near-identity matrix).
