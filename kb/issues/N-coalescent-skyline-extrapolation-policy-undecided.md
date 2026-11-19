# Skyline constant extrapolation policy is unapproved

> [!IMPORTANT]
> **Decision required.** Constant boundary extrapolation is the current runtime behavior, but it is not an approved scientific policy. The options are under "Design axis: out-of-domain evaluation": reject out-of-domain evaluation, keep constant extrapolation and record it in output metadata, or extrapolate with the terminal slope.

The skyline holds $T_c(t)$ as a step function that is constant on each equal-width segment (see [kb/decisions/coalescent-skyline-convex-log-tc.md](../decisions/coalescent-skyline-convex-log-tc.md)): the `PiecewiseConstantFn` schedule [packages/treetime/src/coalescent/skyline.rs#L62-L65](../../packages/treetime/src/coalescent/skyline.rs#L62-L65) stores only the interior segment boundaries, so `fn PiecewiseConstantFn.eval()` returns the value of the first or last segment outside the span, and $T_c$ is extrapolated as a constant. The coalescent prior and the merger rates evaluate this schedule. No output metadata records that values outside the span are extrapolated.

The runtime behavior is therefore explicit. What remains unresolved is whether constant boundary extrapolation is the approved scientific policy and, if retained, how outputs record that values beyond observed lineage events are extrapolated rather than inferred.

## Design axis: out-of-domain evaluation

- O1. Reject evaluation outside the inferred skyline domain with a typed error. This prevents extrapolated values from being mistaken for inferred population history.
- O2. Use constant boundary extrapolation and record the assumption in output metadata.
- O3. Extrapolate with the terminal slope. This can create non-positive or rapidly diverging $T_c$ and requires explicit bounds.

## Recommendation

Use O1 for inference and require an explicit, separately approved extrapolation policy for presentation or downstream prediction.

## Proposed implementation plan

The locations are current: the `PiecewiseConstantFn` schedule at `skyline.rs` L61-L64, and its readers `fn CoalescentPrior::new()`, `fn CoalescentPrior.node_contribution()` and `fn LineageExpectation.merger_rate_schedule()` in `coalescent/mean_field.rs`. The plan depends on the open decision; record the approved option in `kb/decisions/coalescent-skyline-convex-log-tc.md` or a new decision file.

- **O1**: give the schedule its span and return an error for $t$ outside it, and propagate the error through the readers of the schedule. First list every caller that evaluates $T_c$ outside the span (the survival tails of the lineage expectation, see `M-coalescent-lineage-grid-extends-past-the-latest-sample-with-undated-tips.md`), because each must then clip to the span or fail. Tests: evaluation at a time before the first and after the last boundary returns the exact error (`assert_err!`); evaluation inside returns the segment value (oracle: the fitted segment values)
- **O2**: keep the clamping; add a field to the coalescent output (`--output-coalescent`) that states the inferred span and that values outside it are constant extrapolations; run `just gen` for the schemas. Tests: a snapshot of the coalescent output on a small tree that shows the span field (oracle: the segment boundaries of the fit)
- **O3**: extrapolate $\log T_c$ with the slope of the terminal segment pair, bounded by a limit that the decision sets. Tests: values at times outside the span follow the line through the two terminal segment midpoints (oracle: manual calculation)
- Validation for every option: `just check`, `just test-rs`, and `just generated-check` when the output schema changes
