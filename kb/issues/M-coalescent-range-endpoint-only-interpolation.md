# Coalescent weight applied to Range distributions evaluates only at endpoints

`distribution_multiply_by_fn` handles `Distribution::Range` by calling `distribution.t()` and `distribution.y()`, which for a Range return exactly 2 elements (the boundary coordinates and their uniform amplitude). The nonlinear coalescent cost is then evaluated at only those two grid points, and the result is reconstructed as a `Distribution::Function` via linear interpolation.

For any nonlinear weight function -- the coalescent prior's `node_contribution` is nonlinear when $\bar k > 1$ or $T_c$ varies -- this linearly interpolates the output density rather than evaluating the cost at interior coordinates. With $k = 2$, $T_c = 2$, and $H(t) = 2.5 - 0.25t$ on $[0, 10]$, the normalized midpoint should be approximately $0.2865$; endpoint interpolation yields approximately $0.5410$.

## Affected code

- [`packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs#L22-L26`](../../packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs#L22-L26): Range/Point branch extracts 2-element grid via `distribution.t()` and `distribution.y()`
- [`packages/treetime/src/coalescent/mean_field.rs`](../../packages/treetime/src/coalescent/mean_field.rs): `CoalescentPrior::node_contribution` is nonlinear in the time through the integrals of $\bar k/T_c$ and $\ln T_c$ (the `CoalescentModel` contributions this issue first described were removed in phase 2; the endpoint-only evaluation applies to the mean-field prior alike)

## Impact

Uncertain-date nodes whose time constraint is a `Distribution::Range` receive a coalescent factor with incorrect interior values. The error exceeds the project numerical contract ($> 1\mathrm{e}{-6}$) for nontrivial tree configurations and affects inferred node time distributions.

## Constraint

A fix must define a concrete interior grid for nonlinear range weighting without reviving the deleted fine-grid multiplication failure (the original grid-explosion issue). The grid choice is an implementation detail only if it reproduces the v0 oracle within the numerical contract.

## Related issues

- [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md): approved uniform-grid resolution contract for distribution arithmetic
