# Leaf variance offset convention diverges from v0 on terminal edges

## Problem

v1's `EdgeCostFn::evaluate` at `packages/treetime/src/reroot/cost_function.rs:34` computes leaf-edge variance as:

```rust
self.branch_variance * (1.0 - x) + self.variance_offset_leaf
```

With the default `VarianceModel` (`variance_factor=0, variance_offset=0, variance_offset_leaf=1`), `branch_variance = 0` and the effective leaf variance is `0*(1-x) + 1 = 1` for all split fractions `x`. The leaf measurement noise is constant, not scaled by the split position.

v0's `_optimal_root_along_branch` at `packages/legacy/treetime/treetime/treeregression.py:389` passes `var*x` and `var*(1-x)` to `propagate_averages`, where `var = self.branch_variance(n)`. For v0's `min_dev` default (`branch_variance = lambda x: 1.0 if x.is_terminal() else 0.0`), the full `1.0` is scaled by the split fraction: child-side variance is `1.0*(1-x)`, parent-side is `1.0*x`.

## Impact

Affects root positions on **terminal edges only** (2-tip trees, or optima near tips). For interior optima on larger trees -- the normal min-dev case -- both sides have equal unit weights and the objectives coincide.

v1's treatment (fixed tip measurement noise, not split-proportional) is more principled: tip noise represents sequencing error, which does not depend on where on the branch the root sits. v0's treatment (scale everything including the terminal constant by the split fraction) is simpler but conflates branch-length-proportional variance with the constant tip offset.

## Decision needed

- **Match v0**: scale `variance_offset_leaf` by `(1-x)` on the child side and by `x` on the parent side. Simplest parity path.
- **Keep v1**: document as an intentional change in `kb/decisions/`. Scientifically defensible (tip noise is measurement error, not evolutionary distance).

## Code references

- v1 cost function: `packages/treetime/src/reroot/cost_function.rs:30-40`
- v1 variance model: `packages/treetime/src/reroot/variance.rs:11-24`
- v0 split: `packages/legacy/treetime/treetime/treeregression.py:389-392`
- v0 default branch_variance: `packages/legacy/treetime/treetime/clock_tree.py:287`
