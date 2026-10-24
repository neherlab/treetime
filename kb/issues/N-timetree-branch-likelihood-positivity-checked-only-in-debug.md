# Branch length likelihood checks branch rate positivity only in debug builds

`fn compute_branch_length_distribution()` takes one `branch_rate`, the clock rate times the edge's gamma, and asserts its precondition with `debug_assert!(branch_rate > 0.0, ...)` [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L26](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L26). Release builds strip the check. `fn branch_length_grid_to_time_range()` divides the branch-length grid by `branch_rate` [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L99-L104](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L99-L104), so a zero, negative, or NaN value produces an infinite, reversed, or NaN time range without an error.

The only caller, `fn compute_branch_distributions_marginal_mode()`, passes `clock_rate * gamma` (`packages/treetime/src/timetree/inference/runner.rs`). Every producer of these factors keeps them positive:

- **Estimated rate**: when the root-to-tip regression gives a non-positive rate, `timetree` starts from the rate of the Gaussian timetree instead (`fn ClockReroot::into_clock_fit_or()` in `packages/treetime/src/clock/clock_regression.rs` with `fn gaussian_start_rate()` in `packages/treetime/src/timetree/init.rs`, since f6ca696c1), and stops with an error when that rate is undefined (`fn GaussianTree::rate()` in `packages/treetime/src/clock/gaussian.rs`). `fn ClockModel::from_regression_allow_negative()` keeps a non-positive rate only for the `clock` command; `--allow-negative-rate` of `timetree` only lets the root search accept roots with a negative rate
- **Fixed rate**: `fn ClockModel::with_fixed_rate()` rejects a non-positive `--clock-rate` and the non-positive fallback rate (`packages/treetime/src/clock/clock_model.rs`)
- **Marginal clock-rate step**: requires a positive, finite reference rate and multiplies it by an exponential factor (`fn step_marginal_clock_rate()` in `packages/treetime/src/timetree/optimization/marginal_clock_rate.rs`)
- **Gammas**: `fn unit_gammas()` sets every gamma to `1.0`, and `fn apply_relaxed_clock()` floors gamma at `0.1` or falls back to `1.0` (`packages/treetime/src/timetree/optimization/relaxed_clock.rs`)
- **Rate susceptibility of the confidence intervals**: the lower rate is at least 0.1 times the current rate (`fn compute_rate_susceptibility()` in `packages/treetime/src/timetree/confidence.rs`)

The defect is therefore latent: the function's contract depends on these producers, and a new gamma or clock-rate source would fail silently in release builds.

## Impact

- Silent inf/NaN propagation into branch-length distributions in release builds if a non-positive or NaN branch rate reaches the function

## Fix

Replace the `debug_assert!` with a check that runs in all builds and returns a contextual error. `compute_branch_length_distribution()` already returns `Result`, so the error propagates without an API change.

## Validation

- Unit tests for zero, negative, and NaN `branch_rate` that expect an error in both debug and release builds

## Related issues

- [N-error-suppression-unwrap-or-defaults.md](N-error-suppression-unwrap-or-defaults.md): other `debug_assert!` checks that release builds strip
- [N-numerical-stability-magic-constants.md](N-numerical-stability-magic-constants.md)
