# Coalescent weight applied to Range distributions evaluates only at endpoints

`distribution_apply_neg_log_weight` handles `Distribution::Range` by calling `distribution.t()` and `distribution.y()`, which for a Range return exactly 2 elements (the boundary coordinates and their uniform amplitude). The nonlinear coalescent cost is then evaluated at only those two grid points, and the result is reconstructed as a `Distribution::Function` via linear interpolation.

For any nonlinear weight function -- the coalescent `internal_contribution` and `root_contribution` are nonlinear when $k > 1$ or $T_c$ varies -- this linearly interpolates the output density rather than evaluating the cost at interior coordinates. With $k = 2$, $T_c = 2$, and $H(t) = 2.5 - 0.25t$ on $[0, 10]$, the normalized midpoint should be approximately $0.2865$; endpoint interpolation yields approximately $0.5410$.

## Affected code

- [`packages/treetime-distribution/src/distribution_ops/log_cost.rs#L33-L36`](../../packages/treetime-distribution/src/distribution_ops/log_cost.rs#L33-L36): Range/Point branch extracts 2-element grid via `distribution.t()` and `distribution.y()`
- [`packages/treetime/src/coalescent/coalescent.rs#L65-L73`](../../packages/treetime/src/coalescent/coalescent.rs#L65-L73): `internal_contribution` and `root_contribution` are nonlinear in $H(t)$ and $\ln(\lambda(t))$

## Impact

Uncertain-date nodes whose time constraint is a `Distribution::Range` receive a coalescent factor with incorrect interior values. The error exceeds the project numerical contract ($> 1\mathrm{e}{-6}$) for nontrivial tree configurations and affects inferred node time distributions.

## Constraint

A fix must define a concrete interior grid for nonlinear range weighting without reviving the deleted fine-grid multiplication failure (the original grid-explosion issue). The grid choice is an implementation detail only if it reproduces the v0 oracle within the numerical contract.

## Related issues

- [M-distribution-product-grid-resolution-diverges-from-v0.md](M-distribution-product-grid-resolution-diverges-from-v0.md): grid resolution divergence in distribution products (separate concern, same domain)
- [M-timetree-backward-pass-plain-space-underflow.md](M-timetree-backward-pass-plain-space-underflow.md): underflow in the same backward pass product path (separate mechanism)
