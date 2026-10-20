# Skyline constant extrapolation policy is unapproved

> [!IMPORTANT]
> **Decision required.** Constant boundary extrapolation is the current runtime behavior, but it is not an approved scientific policy. The options are under "Design axis: out-of-domain evaluation": reject out-of-domain evaluation, keep constant extrapolation and record it in output metadata, or extrapolate with the terminal slope.

`build_tc_distribution()` [packages/treetime/src/coalescent/skyline.rs#L298-L314](../../packages/treetime/src/coalescent/skyline.rs#L298-L314) builds $T_c(t)$ as a step function (a `Distribution::Formula` that is constant on each equal-width segment, see [kb/decisions/coalescent-skyline-convex-log-tc.md](../decisions/coalescent-skyline-convex-log-tc.md)). Outside the segment span, `segment_index()` clamps to the first or last segment, so $T_c$ is extrapolated as a constant. The `PiecewiseConstantFn` schedule of the skyline [packages/treetime/src/coalescent/skyline.rs#L63-L66](../../packages/treetime/src/coalescent/skyline.rs#L63-L66) does the same, because it stores only the interior segment boundaries. No output metadata records that values outside the span are extrapolated.

The runtime behavior is therefore explicit. What remains unresolved is whether constant boundary extrapolation is the approved scientific policy and, if retained, how outputs record that values beyond observed lineage events are extrapolated rather than inferred.

## Design axis: out-of-domain evaluation

- O1. Reject evaluation outside the inferred skyline domain with a typed error. This prevents extrapolated values from being mistaken for inferred population history.
- O2. Use constant boundary extrapolation and record the assumption in output metadata.
- O3. Extrapolate with the terminal slope. This can create non-positive or rapidly diverging $T_c$ and requires explicit bounds.

## Recommendation

Use O1 for inference and require an explicit, separately approved extrapolation policy for presentation or downstream prediction.
