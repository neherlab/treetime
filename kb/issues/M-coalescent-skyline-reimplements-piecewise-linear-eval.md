# Skyline build_tc_distribution reimplements PiecewiseLinearFn::eval

`build_tc_distribution` in `coalescent/skyline.rs` contains a closure that performs boundary clamping, `partition_point` lookup, and linear interpolation -- the same algorithm as `PiecewiseLinearFn::eval()` in `treetime-grid`. The skyline module already imports from `treetime-grid` but does not use `PiecewiseLinearFn` for this evaluation.

## Locations

- `packages/treetime/src/coalescent/skyline.rs:230-252:` inline closure reimplementing piecewise linear evaluation over `time_grid` and `tc_values`
- `packages/treetime-grid/src/piecewise_linear_fn.rs:75-95:` `fn eval()` -- identical algorithm (boundary clamp, `partition_point`, linear interpolation)

Both use the same steps: check boundary conditions, find interval via `partition_point`, compute interpolation fraction `alpha = (t - t0) / (t1 - t0)`, return `v0 + alpha * (v1 - v0)`.

## Impact

If the interpolation algorithm is corrected or extended (e.g., extrapolation behavior), the inline version won't receive the fix.

## Action

Build a `PiecewiseLinearFn` from `(time_grid, tc_values)` and wrap its `.eval()` method in the `DistributionFormula` closure.
