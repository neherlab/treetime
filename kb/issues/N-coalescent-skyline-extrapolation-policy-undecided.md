# Skyline constant extrapolation policy is unapproved

> **Still open; implementation changed.** `build_tc_distribution()` now emits a
> piecewise-*constant* $T_c(t)$ (a step function via `Distribution::Formula`), not
> a piecewise-linear `GridFn`. Out-of-domain times still clamp to the first/last
> segment (constant extrapolation), so the policy question below remains. See
> [decisions/coalescent-skyline-convex-inverse-tc.md](../decisions/coalescent-skyline-convex-inverse-tc.md).

`build_tc_distribution()` constructs a piecewise-constant $T_c(t)$ distribution and clamps to the first/last segment value outside the segment grid.

The runtime behavior is therefore explicit. What remains unresolved is whether constant boundary extrapolation is the approved scientific policy and, if retained, how outputs record that values beyond observed lineage events are extrapolated rather than inferred.

## Design axis: out-of-domain evaluation

- O1. Reject evaluation outside the inferred skyline domain with a typed error. This prevents extrapolated values from being mistaken for inferred population history.
- O2. Use constant boundary extrapolation and record the assumption in output metadata.
- O3. Extrapolate with the terminal slope. This can create non-positive or rapidly diverging $T_c$ and requires explicit bounds.

## Recommendation

Use O1 for inference and require an explicit, separately approved extrapolation policy for presentation or downstream prediction.

## Related issues

- [N-coalescent-skyline-grid-validation-incomplete.md](N-coalescent-skyline-grid-validation-incomplete.md)
