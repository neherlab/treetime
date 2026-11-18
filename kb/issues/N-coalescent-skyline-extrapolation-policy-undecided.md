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

## Proposed implementation plan

The locations are current: `fn build_tc_distribution()` at `skyline.rs` L298, `fn segment_index()` at L310 (clamps to the first or last segment), and the `PiecewiseConstantFn` schedule at L63-L66. The plan depends on the open decision; record the approved option in `kb/decisions/coalescent-skyline-convex-log-tc.md` or a new decision file.

- **O1**: make `fn segment_index()` return an error for $t$ outside the segment span, and propagate it through the $T_c$ formula and the `PiecewiseConstantFn` schedule. First list every caller that evaluates $T_c$ outside the span (the survival tails of the lineage expectation, see `M-coalescent-lineage-grid-extends-past-the-latest-sample-with-undated-tips.md`), because each must then clip to the span or fail. Tests: evaluation at a time before the first and after the last boundary returns the exact error (`assert_err!`); evaluation inside returns the segment value (oracle: the fitted segment values)
- **O2**: keep the clamping; add a field to the coalescent output (`--output-coalescent`) that states the inferred span and that values outside it are constant extrapolations; run `just gen` for the schemas. Tests: a snapshot of the coalescent output on a small tree that shows the span field (oracle: the segment boundaries of the fit)
- **O3**: extrapolate $\log T_c$ with the slope of the terminal segment pair, bounded by a limit that the decision sets. Tests: values at times outside the span follow the line through the two terminal segment midpoints (oracle: manual calculation)
- Validation for every option: `just check`, `just test-rs`, and `just generated-check` when the output schema changes
