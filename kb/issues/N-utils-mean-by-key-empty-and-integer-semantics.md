# `mean_by_key` returns zero for empty input and truncates integer means

## Summary

`MeanByKey::mean_by_key` in `packages/treetime-utils/src/iterator/mean_by_key.rs` has two behaviors that the unit tests in `packages/treetime-utils/src/iterator/__tests__/test_mean_by_key.rs` fix in place without stating a reason:

- **Empty input returns zero**: an empty iterator yields `U::zero()`, so "no data" and "mean is zero" are indistinguishable
- **Integer means truncate**: for integer key types the division truncates, so the mean of `1..=1000` is `500` instead of `500.5`
- **Count conversion falls back to zero**: when `U::from_usize(count)` fails, the divisor becomes `U::zero()`, which divides by zero for integer types

## Current callers

The only production callers are in `packages/treetime-validation/src/testing/console/console_metrics.rs`. They use `f64` keys over groups that are non-empty by construction, so none of these behaviors affects current output.

## Open question

Decide the contract before a caller relies on it:

- **Return `Option<U>`**: `None` for empty input, and an error or `None` when the count does not fit `U`
- **Restrict to floating-point keys**: remove integer support, so truncation cannot occur
- **Keep the current contract**: document zero for empty input and integer truncation as intended behavior
