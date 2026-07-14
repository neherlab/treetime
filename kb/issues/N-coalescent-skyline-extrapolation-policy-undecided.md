# Skyline extrapolation policy is implicit

`fn build_tc_distribution()` constructs a piecewise-linear $T_c(t)$ distribution over the time grid [packages/treetime/src/coalescent/skyline.rs#L215-L228](../../packages/treetime/src/coalescent/skyline.rs#L215-L228). Evaluation outside that grid inherits boundary behavior from `PiecewiseLinearFn`, so the scientific assumption beyond observed lineage events is not declared at the skyline boundary.

## Design axis: out-of-domain evaluation

- O1. Reject evaluation outside the inferred skyline domain with a typed error. This prevents extrapolated values from being mistaken for inferred population history.
- O2. Use constant boundary extrapolation and record the assumption in output metadata.
- O3. Extrapolate with the terminal slope. This can create non-positive or rapidly diverging $T_c$ and requires explicit bounds.

## Recommendation

Use O1 for inference and require an explicit, separately approved extrapolation policy for presentation or downstream prediction.

## Related issues

- [N-coalescent-skyline-grid-validation-incomplete.md](N-coalescent-skyline-grid-validation-incomplete.md)
