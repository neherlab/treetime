# Distribution product grid resolution diverges from v0

V0 and v1 preserve different information when multiplying discretized distributions. The difference affects coalescent ordering comparisons and any operation whose operands have non-aligned knots.

## Reference behavior

V0 [`Distribution.multiply()`](../../packages/legacy/treetime/treetime/distribution.py) constructs the sorted union of all input knots, clips that union to the support intersection, and evaluates every operand on the consolidated grid. For products of more than three distributions it can thin the center of a very large union, but peripheral knots remain explicit.

V1 [`multiply_function_function()`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs) discards both input grids and creates a uniform grid with `max(a.len(), b.len())` points over the support intersection. `multiply_formula_function()` uses the current function length, while Formula-only operations use a fixed `FORMULA_GRID_SIZE`. Sequential multiplication can therefore shift or remove knots on every step, and the result depends on multiplication order.

## Unresolved accuracy contract

There is no stated error bound connecting the chosen uniform point count to interpolation error, posterior peak displacement, or integrated probability. The `MAX_GRID_POINTS` cap proposed in commit `542ac860c7cfa4bab6764aee1d1b3810a09eb54f` likewise has no derived approximation-error bound. A cap could prevent pathological growth while silently discarding sharp features.

No implementation ticket is ready until the project chooses one of these contracts and its validation criterion:

- preserve v0's union-of-knots behavior, including a parity-compatible thinning rule;
- use a derived spacing/error criterion for uniform resampling; or
- approve a bounded approximation with empirically justified error limits.

## Related issues

- [H-timetree-coalescent-grid-explosion.md](H-timetree-coalescent-grid-explosion.md)
- [M-timetree-coalescent-multiplication-ordering-diverges-from-v0.md](M-timetree-coalescent-multiplication-ordering-diverges-from-v0.md)
- [M-distribution-mixed-support-operations-lose-exact-intersection-boundaries.md](M-distribution-mixed-support-operations-lose-exact-intersection-boundaries.md)
