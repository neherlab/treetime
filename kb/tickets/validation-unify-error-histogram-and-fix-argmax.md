# Unify error histogram functions and replace argmax with argmax_first

Two quality issues in validation metrics: duplicated histogram logic and non-deterministic argmax.

## Current state

`histograms.rs:75-133` (`compute_error_histogram`) and `:135-177` (`compute_error_histogram_signed`) share 80% logic, differing in filter predicate and log-scale option. `peak_metrics.rs:27-28` uses `ndarray_stats::argmax()` instead of project's `argmax_first()`.

## Target state

A single `compute_error_histogram(errors, config)` function with a config enum or parameter controlling signed/unsigned mode. `argmax()` replaced with `argmax_first()`.

## Implementation

1. Define `enum HistogramMode { Unsigned, Signed }` or equivalent parameter
2. Merge `compute_error_histogram` and `compute_error_histogram_signed` into one function parameterized by mode
3. Mode controls: filter predicate (`x > 0.0` vs `x.is_finite()`) and value transform (`log10` vs identity)
4. Update all callers to pass the appropriate mode
5. Replace `argmax()` with `argmax_first()` in `peak_metrics.rs:27-28`
6. Verify validation test output unchanged

## Related issues

Source: `kb/issues/N-validation-error-histogram-and-argmax-quality.md`
