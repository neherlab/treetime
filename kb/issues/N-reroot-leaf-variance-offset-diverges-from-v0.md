# Leaf variance offset convention diverges from v0 on terminal edges

> [!IMPORTANT]
> **Decision required.** The project rule makes exact v0 parity the default, and a deviation counts as a defect until it is approved. That rule points to "Match v0". This issue also asks for an empirical comparison of both conventions before either is chosen, because v1's fixed tip-noise interpretation is scientifically defensible. The v1 behavior is not an approved divergence. Decide between "Match v0" and "Keep v1" (with a `kb/decisions/` entry), or confirm that parity applies without the comparison.

## Problem

v1's `pub(crate) fn EdgeCostFn::evaluate()` [packages/treetime/src/reroot/cost_function.rs#L27-L45](../../packages/treetime/src/reroot/cost_function.rs#L27-L45) computes leaf-edge variance as:

```rust
self.branch_variance * (1.0 - x) + self.variance_offset_leaf
```

With the default `VarianceModel` (`variance_factor=0, variance_offset=0, variance_offset_leaf=1`), `branch_variance = 0` and the effective leaf variance is `0*(1-x) + 1 = 1` for all split fractions `x`. The leaf measurement noise is constant, not scaled by the split position.

The clock reroot search uses the same convention: `pub(crate) fn BranchPointCostFunction::evaluate_clock_set()` [packages/treetime/src/clock/find_best_root/cost_function.rs#L74-L80](../../packages/treetime/src/clock/find_best_root/cost_function.rs#L74-L80) adds `self.options.variance_offset_leaf` unscaled to the child side.

v0's `def TreeRegression._optimal_root_along_branch()` [packages/legacy/treetime/treetime/treeregression.py#L385](../../packages/legacy/treetime/treetime/treeregression.py#L385) passes `var*x` and `var*(1-x)` to `propagate_averages`, where `var = self.branch_variance(n)`. For v0's `min_dev` default [packages/legacy/treetime/treetime/clock_tree.py#L287](../../packages/legacy/treetime/treetime/clock_tree.py#L287), the full `1.0` is scaled by the split fraction: child-side variance is `1.0*(1-x)`, parent-side is `1.0*x`.

## Impact

Affects root positions on **terminal edges only** (2-tip trees, or optima near tips). For interior optima on larger trees -- the normal min-dev case -- both sides have equal unit weights and the objectives coincide.

V1's fixed tip-noise interpretation is scientifically plausible, but it is not approved. Exact v0 parity remains the required default.

## Options

- **Match v0**: scale `variance_offset_leaf` by `(1-x)` on the child side and by `x` on the parent side. Simplest parity path.
- **Keep v1**: document as an intentional change in `kb/decisions/`. Scientifically defensible (tip noise is measurement error, not evolutionary distance).

## Recommendation

Match v0 by splitting the complete terminal variance between both sides. Retaining v1 requires explicit consent and a decision entry; the issue itself does not grant that consent.

## Fix (Match v0)

- Compute the full terminal variance, including the leaf offset, once
- Apply the fractions $x$ and $1-x$ to the two sides of the candidate root split
- Use the same convention in the clock, timetree, and optimize reroot paths: the generic `EdgeCostFn` (optimize) and the clock `BranchPointCostFunction` (clock and timetree)

## Validation

- Side variances at $x=0$, $x=0.5$, and $x=1$
- Golden masters against v0 for a two-tip tree and for an optimum near a tip
- The clock and optimize paths give the same result for the same objective and inputs

## Code references

- `pub(crate) fn EdgeCostFn::evaluate()` [packages/treetime/src/reroot/cost_function.rs#L27](../../packages/treetime/src/reroot/cost_function.rs#L27)
- `pub(crate) fn BranchPointCostFunction::evaluate_clock_set()` [packages/treetime/src/clock/find_best_root/cost_function.rs#L74](../../packages/treetime/src/clock/find_best_root/cost_function.rs#L74)
- `struct VarianceModel` [packages/treetime/src/reroot/variance.rs#L5](../../packages/treetime/src/reroot/variance.rs#L5)
- `def TreeRegression._optimal_root_along_branch()` [packages/legacy/treetime/treetime/treeregression.py#L385](../../packages/legacy/treetime/treetime/treeregression.py#L385)
- v0 `min_dev` default [packages/legacy/treetime/treetime/clock_tree.py#L287](../../packages/legacy/treetime/treetime/clock_tree.py#L287)

## Related KB items

- [kb/proposals/reroot-generic-scoring-architecture.md](../proposals/reroot-generic-scoring-architecture.md)
- [kb/proposals/optimize-reroot-support.md](../proposals/optimize-reroot-support.md)
- [N-reroot-missing-min-dev-end-to-end-oracle.md](N-reroot-missing-min-dev-end-to-end-oracle.md)
- [N-reroot-clock-search-duplicates-generic-module.md](N-reroot-clock-search-duplicates-generic-module.md)
