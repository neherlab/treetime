# Coalescent weight applied to Range distributions evaluates only at endpoints

`distribution_multiply_by_fn` handles `Distribution::Range` by calling `distribution.t()` and `distribution.y()`, which for a Range return exactly 2 elements (the boundary coordinates and their uniform amplitude). The nonlinear coalescent cost is then evaluated at only those two grid points, and the result is reconstructed as a `Distribution::Function` via linear interpolation.

For any nonlinear weight function -- the coalescent prior's `node_contribution` is nonlinear when $\bar k > 1$ or $T_c$ varies -- this linearly interpolates the output density rather than evaluating the cost at interior coordinates. With $k = 2$, $T_c = 2$, and $H(t) = 2.5 - 0.25t$ on $[0, 10]$, the normalized midpoint should be approximately $0.2865$; endpoint interpolation yields approximately $0.5410$.

## Affected code

- [`packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs#L22-L30`](../../packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs#L22-L30): the Range/Point branch extracts the 2-element grid via `distribution.t()` and `distribution.y()` and evaluates the weight on it
- [`packages/treetime/src/coalescent/mean_field.rs`](../../packages/treetime/src/coalescent/mean_field.rs): `CoalescentPrior::node_contribution` is nonlinear in the time through the integrals of $\bar k/T_c$ and $\ln T_c$

## Impact

No production path reaches the endpoint-only evaluation. The only caller that can pass a Range is the coalescent factor of tips, and `fn apply_leaf_coalescent_factor()` [packages/treetime/src/timetree/inference/backward_pass.rs#L135-L158](../../packages/treetime/src/timetree/inference/backward_pass.rs#L135-L158) first converts a Range date into a `Distribution::Function` with hard edges on a grid at the time resolution of the coalescent prior (`time_resolution()`, at least 2 and at most `GRID_POINTS` points). The coalescent factor of internal nodes, `fn apply_coalescent_prior()` [packages/treetime/src/timetree/inference/backward_pass.rs#L106-L119](../../packages/treetime/src/timetree/inference/backward_pass.rs#L106-L119), receives the product of the child messages.

The defect remains in the API of `distribution_multiply_by_fn()`: a direct call with a Range gives a result whose interior values are wrong by more than the project numerical contract ($> 1\mathrm{e}{-6}$) for a nonlinear weight. The remaining work is to reject Range input in `distribution_multiply_by_fn()` with an error, or to grid it there as `apply_leaf_coalescent_factor()` does, so the gridding has one owner.

## Constraint

A fix must define a concrete interior grid for nonlinear range weighting without reviving the deleted fine-grid multiplication failure (the original grid-explosion issue). The grid choice is an implementation detail only if it reproduces the v0 oracle within the numerical contract.

## Related issues

- [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md): approved uniform-grid resolution contract for distribution arithmetic
