# `mean_by_key` returns zero for empty input and truncates integer means

## Summary

`MeanByKey::mean_by_key` in `packages/treetime-utils/src/iterator/mean_by_key.rs` has three behaviors that the unit tests in `packages/treetime-utils/src/iterator/__tests__/test_mean_by_key.rs` fix in place without stating a reason:

- **Empty input returns zero**: an empty iterator yields `U::zero()`, so "no data" and "mean is zero" are indistinguishable
- **Integer means truncate**: for integer key types the division truncates, so the mean of `1..=1000` is `500` instead of `500.5`
- **Count conversion falls back to zero**: when `U::from_usize(count)` fails, the divisor becomes `U::zero()`, which divides by zero for integer types

## Current callers

The only production callers are in `packages/treetime-validation/src/testing/console/console_metrics.rs`. They use `f64` keys over groups that are non-empty by construction, so none of these behaviors affects current output.

## Open question

> [!IMPORTANT]
> **Decision required.** The contract of `mean_by_key()` is undefined for empty input, integer keys, and counts that do not fit `U`. No current caller depends on it, so the choice is free. Options:
>
> - **Return `Option<U>`**: `None` for empty input, and an error or `None` when the count does not fit `U`
> - **Restrict to floating-point keys**: remove integer support, so truncation cannot occur
> - **Keep the current contract**: document zero for empty input and integer truncation as intended behavior

## Proposed implementation plan

The plan depends on the open decision; each option changes `MeanByKey::mean_by_key()` in `packages/treetime-utils/src/iterator/mean_by_key.rs`, its tests in `packages/treetime-utils/src/iterator/__tests__/test_mean_by_key.rs`, and the callers in `packages/treetime-validation/src/testing/console/console_metrics.rs`.

- **Return `Option<U>`**: return `None` for empty input and for a count that `U::from_usize()` cannot represent; callers unwrap with a message where the group is non-empty by construction. Tests: empty -> `None`, `[1.0, 2.0]` -> `Some(1.5)`; oracle: hand calculation
- **Restrict to floating-point keys**: bound `U: Float` (`num_traits`), convert the count with `NumCast` and propagate a failure; tests over `f64` only, `1..=1000` as `f64` -> `500.5` (oracle: (1 + 1000) / 2)
- **Keep the current contract**: replace `unwrap_or_else(U::zero)` with an error or a documented panic-free path, and keep the tests

Validate with `just test-rs -E 'package(treetime-utils)'` and `-E 'package(treetime-validation)'`.
