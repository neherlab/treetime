# Branch length likelihood checks clock rate and gamma positivity only in debug builds

`fn compute_branch_length_distribution()` asserts its preconditions with `debug_assert!(clock_rate > 0.0, ...)` and `debug_assert!(gamma > 0.0)` [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L26-L27](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L26-L27). Release builds strip both checks. The function divides the branch-length grid by `clock_rate * gamma` [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L103-L109](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L103-L109), so a zero, negative, or NaN value produces an infinite, reversed, or NaN time range without an error.

The current producers keep both values positive:

- The clock model rejects a non-positive estimated or specified clock rate [packages/treetime/src/clock/clock_model.rs#L45](../../packages/treetime/src/clock/clock_model.rs#L45), [#L62](../../packages/treetime/src/clock/clock_model.rs#L62), [#L89](../../packages/treetime/src/clock/clock_model.rs#L89)
- `fn unit_gammas()` sets every gamma to `1.0`, and `fn apply_relaxed_clock()` clamps gamma to at least `0.1` or falls back to `1.0` [packages/treetime/src/timetree/optimization/relaxed_clock.rs#L65-L85](../../packages/treetime/src/timetree/optimization/relaxed_clock.rs#L65-L85)
- The confidence-interval runs scale gamma by a positive rate ratio [packages/treetime/src/timetree/confidence.rs#L32-L40](../../packages/treetime/src/timetree/confidence.rs#L32-L40)

The defect is therefore latent: the function's contract depends on these producers, and a new gamma or clock-rate source would fail silently in release builds.

## Impact

- Silent inf/NaN propagation into branch-length distributions in release builds if a non-positive or NaN clock rate or gamma reaches the function

## Fix

Replace both `debug_assert!` calls with checks that run in all builds and return a contextual error. `compute_branch_length_distribution()` already returns `Result`, so the error propagates without an API change.

## Validation

- Unit tests for zero, negative, and NaN `clock_rate` and `gamma` that expect an error in both debug and release builds

## Related issues

- [N-error-suppression-unwrap-or-defaults.md](N-error-suppression-unwrap-or-defaults.md): other `debug_assert!` checks that release builds strip
- [N-numerical-stability-magic-constants.md](N-numerical-stability-magic-constants.md)
