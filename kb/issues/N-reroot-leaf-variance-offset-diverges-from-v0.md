# Leaf variance offset convention diverges from v0 on terminal edges

> [!IMPORTANT]
> **Research and discussion required.** The v1 behavior is not an approved divergence. Compare both conventions empirically and obtain an explicit decision before changing or retaining it.

## Problem

v1's `pub fn EdgeCostFn::evaluate()` [packages/treetime/src/reroot/cost_function.rs#L29](../../packages/treetime/src/reroot/cost_function.rs#L29) computes leaf-edge variance as:

```rust
self.branch_variance * (1.0 - x) + self.variance_offset_leaf
```

With the default `VarianceModel` (`variance_factor=0, variance_offset=0, variance_offset_leaf=1`), `branch_variance = 0` and the effective leaf variance is `0*(1-x) + 1 = 1` for all split fractions `x`. The leaf measurement noise is constant, not scaled by the split position.

v0's `def TreeRegression._optimal_root_along_branch()` [packages/legacy/treetime/treetime/treeregression.py#L385](../../packages/legacy/treetime/treetime/treeregression.py#L385) passes `var*x` and `var*(1-x)` to `propagate_averages`, where `var = self.branch_variance(n)`. For v0's `min_dev` default [packages/legacy/treetime/treetime/clock_tree.py#L287](../../packages/legacy/treetime/treetime/clock_tree.py#L287), the full `1.0` is scaled by the split fraction: child-side variance is `1.0*(1-x)`, parent-side is `1.0*x`.

## Impact

Affects root positions on **terminal edges only** (2-tip trees, or optima near tips). For interior optima on larger trees -- the normal min-dev case -- both sides have equal unit weights and the objectives coincide.

V1's fixed tip-noise interpretation is scientifically plausible, but it is not approved. Exact v0 parity remains the required default.

## Decision needed

- **Match v0**: scale `variance_offset_leaf` by `(1-x)` on the child side and by `x` on the parent side. Simplest parity path.
- **Keep v1**: document as an intentional change in `kb/decisions/`. Scientifically defensible (tip noise is measurement error, not evolutionary distance).

## Recommendation

Match v0 by splitting the complete terminal variance between both sides. Retaining v1 requires explicit consent and a decision entry; the issue itself does not grant that consent.

## Code references

- `pub fn EdgeCostFn::evaluate()` [packages/treetime/src/reroot/cost_function.rs#L29](../../packages/treetime/src/reroot/cost_function.rs#L29)
- `struct VarianceModel` [packages/treetime/src/reroot/variance.rs#L12](../../packages/treetime/src/reroot/variance.rs#L12)
- `def TreeRegression._optimal_root_along_branch()` [packages/legacy/treetime/treetime/treeregression.py#L385](../../packages/legacy/treetime/treetime/treeregression.py#L385)
- v0 `min_dev` default [packages/legacy/treetime/treetime/clock_tree.py#L287](../../packages/legacy/treetime/treetime/clock_tree.py#L287)

## Related KB items

- [kb/proposals/reroot-generic-scoring-architecture.md](../proposals/reroot-generic-scoring-architecture.md)
- [kb/proposals/optimize-reroot-support.md](../proposals/optimize-reroot-support.md)
