# Replace inline interpolation in skyline with PiecewiseLinearFn

`build_tc_distribution` reimplements `PiecewiseLinearFn::eval()` as an inline closure instead of using the existing struct from `treetime-grid`.

## Current state

`coalescent/skyline.rs:230-252` contains a closure that performs boundary clamping, `partition_point` lookup, and linear interpolation over `time_grid` and `tc_values`. `treetime-grid::PiecewiseLinearFn::eval()` implements the identical algorithm.

## Target state

`build_tc_distribution` constructs a `PiecewiseLinearFn` from `(time_grid, tc_values)` and wraps `.eval()` in the `DistributionFormula` closure.

## Implementation

1. In `build_tc_distribution`, construct `let pwl = PiecewiseLinearFn::new(time_grid.clone(), tc_values.clone())?;`
2. Replace the inline interpolation closure with `move |t| Ok(pwl.eval(t))`
3. Verify that boundary behavior matches (both clamp to first/last value)
4. Run coalescent tests to confirm identical results

## Related issues

Source: `kb/issues/M-coalescent-skyline-reimplements-piecewise-linear-eval.md`
