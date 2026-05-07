# ClockSet::chisq() numerically unstable for near-zero determinant

`ClockSet::chisq()` divides by `determinant()` without a minimum threshold guard. When `det` is positive but near zero, the chi-squared value can become extreme (large negative), winning the `< best_chisq` comparison in `find_best_root()` and selecting a near-degenerate root position.

## Impact

In practice, with `ClockParams::default()` (zero `variance_factor`), the determinant at internal nodes is `n^2 * var(dates)`, which is well-conditioned for any dataset with date variation. The risk is theoretical.

The `force_positive_rate: false` path in the pre-filter step does not meaningfully increase exposure: `force_positive=true` already accepts any `det > 0` (no minimum threshold), and the additional nodes that pass with `force_positive=false` are those with negative rate, not near-zero determinant.

## v1 Location

`packages/treetime/src/commands/clock/clock_set.rs:118-125` (`chisq()`)

`packages/treetime/src/commands/clock/find_best_root/find_best_root.rs:137-139` (`has_positive_clock_rate` checks `det > 0.0` without minimum)

## Related issues

- Source: [N-clock-chisq-near-zero-determinant.md](../issues/N-clock-chisq-near-zero-determinant.md) -- delete after full resolution
