# Optimizer evaluation functions omit mu scaling factor

The branch-length evaluation functions (`evaluate_dense_contribution_impl`, `evaluate_sparse_contribution_impl`) compute `exp(eigvals * t)` instead of `exp(mu * eigvals * t)`. The GTR `expQt()` method uses `exp(mu * eigvals * t)`. This discrepancy means the optimizer and `expQt()` interpret branch lengths in different scales when `mu != 1`.

## Impact

For JC69 with default `mu = 1.0`: no impact. For non-unit `mu` models: optimized branch lengths are in the wrong scale by a factor of `mu`. The `is_zero_branch_optimal()` derivative-sign criterion is unaffected because `mu > 0` does not change the sign of the weighted eigenvalue average.

All evaluation code (dense, sparse, unified, zero-branch) is internally consistent: all omit `mu`. The error manifests only when optimized branch lengths are used with `expQt()` or compared across different `mu` values.

## Proposed solution

Multiply eigenvalues by `mu` in all evaluation functions, or absorb `mu` into the eigenvalues at GTR construction time.
