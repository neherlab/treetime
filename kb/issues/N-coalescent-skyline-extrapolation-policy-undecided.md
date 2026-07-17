# Skyline constant extrapolation policy is unapproved

`fn build_tc_distribution()` constructs a piecewise-linear $T_c(t)$ distribution over the time grid [packages/treetime/src/coalescent/skyline.rs#L208-L216](../../packages/treetime/src/coalescent/skyline.rs#L208-L216). The underlying `GridFn` returns the first or last grid value outside that domain [packages/treetime-grid/src/grid_fn.rs#L267-L335](../../packages/treetime-grid/src/grid_fn.rs#L267-L335), and a production-level regression test verifies this constant extrapolation [packages/treetime/src/coalescent/__tests__/test_skyline.rs#L173-L200](../../packages/treetime/src/coalescent/__tests__/test_skyline.rs#L173-L200).

The runtime behavior is therefore explicit. What remains unresolved is whether constant boundary extrapolation is the approved scientific policy and, if retained, how outputs record that values beyond observed lineage events are extrapolated rather than inferred.

## Design axis: out-of-domain evaluation

- O1. Reject evaluation outside the inferred skyline domain with a typed error. This prevents extrapolated values from being mistaken for inferred population history.
- O2. Use constant boundary extrapolation and record the assumption in output metadata.
- O3. Extrapolate with the terminal slope. This can create non-positive or rapidly diverging $T_c$ and requires explicit bounds.

## Recommendation

Use O1 for inference and require an explicit, separately approved extrapolation policy for presentation or downstream prediction.

## Related issues

- [N-coalescent-skyline-grid-validation-incomplete.md](N-coalescent-skyline-grid-validation-incomplete.md)
