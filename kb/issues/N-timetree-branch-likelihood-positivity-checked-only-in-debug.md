# Branch length likelihood checks branch rate positivity only in debug builds

`fn compute_branch_length_distribution()` takes one `branch_rate`, the clock rate times the edge's gamma, and asserts its precondition with `debug_assert!(branch_rate > 0.0, ...)` [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L26](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L26). Release builds strip the check. `fn branch_length_grid_to_time_range()` divides the branch-length grid by `branch_rate` [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L99-L104](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L99-L104), so a zero, negative, or NaN value produces an infinite, reversed, or NaN time range without an error.

The only caller, `fn compute_branch_distributions_marginal_mode()`, passes `clock_rate * gamma` [packages/treetime/src/timetree/inference/runner.rs#L280](../../packages/treetime/src/timetree/inference/runner.rs#L280). The current producers keep both factors positive:

- The timetree clock model rejects a non-positive estimated rate in `fn ClockModel::from_regression()` [packages/treetime/src/clock/clock_model.rs#L42-L57](../../packages/treetime/src/clock/clock_model.rs#L42-L57) and a non-positive specified rate in `fn ClockModel::with_fixed_rate()` [packages/treetime/src/clock/clock_model.rs#L85-L100](../../packages/treetime/src/clock/clock_model.rs#L85-L100). `fn ClockModel::from_regression_allow_negative()` keeps a non-positive rate, but only the `clock` command uses it [packages/treetime/src/clock/pipeline.rs#L170](../../packages/treetime/src/clock/pipeline.rs#L170)
- The marginal clock-rate step multiplies the positive reference rate by an exponential factor [packages/treetime/src/timetree/optimization/marginal_clock_rate.rs#L39](../../packages/treetime/src/timetree/optimization/marginal_clock_rate.rs#L39)
- `fn unit_gammas()` sets every gamma to `1.0`, and `fn apply_relaxed_clock()` clamps gamma to at least `0.1` or falls back to `1.0` [packages/treetime/src/timetree/optimization/relaxed_clock.rs#L65-L85](../../packages/treetime/src/timetree/optimization/relaxed_clock.rs#L65-L85)
- The rate-susceptibility runs of the confidence intervals scale gamma by a positive rate ratio [packages/treetime/src/timetree/confidence.rs#L33-L41](../../packages/treetime/src/timetree/confidence.rs#L33-L41)

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
