# Confidence interval output not gated on --confidence

v1 writes `confidence_intervals.tsv` whenever `--time-marginal` is `always` or `only-final`, regardless of whether `--confidence` is also passed. v0 only writes CI output when `--confidence` is explicitly requested.

## v0 behavior

CI file output is gated on the `calc_confidence` flag (wrappers.py:580, wrappers.py:905-911). The `--confidence` flag must be set AND prerequisites (`--covariation` or `--clock-std-dev`) must be satisfied for CIs to appear in output. Without `--confidence`, no CI data is written even when marginal distributions are available from `time_marginal='always'`.

## v1 behavior

The CI-writing condition checks `time_marginal` mode directly:

```rust
if matches!(time_marginal, TimeMarginalMode::OnlyFinal | TimeMarginalMode::Always) || rate_std.is_some() {
```

This means:

- `--time-marginal=always` or `--time-marginal=only-final` writes mutation-stochasticity CIs (from HPD regions of marginal posterior distributions)
- `--confidence` with prerequisites adds rate-uncertainty CIs (from rate susceptibility analysis), combined via quadrature sum
- Both sources are independent and can appear separately

## Rationale

Marginal posterior distributions contain uncertainty information (HPD regions). Suppressing this information behind a separate `--confidence` flag when the user has already requested marginal mode adds a surprising gate. Users who pass `--time-marginal=always` or `--time-marginal=only-final` are requesting marginal distributions, and confidence intervals are a direct product of those distributions.

The `--confidence` flag in v1 controls rate susceptibility analysis (re-running inference at rate +/- sigma), which is computationally expensive and requires additional prerequisites. Separating "write mutation-based CIs" from "also compute rate-based CIs" gives users more granular control.

## Impact

Users running `--time-marginal=always` or `--time-marginal=only-final` get a `confidence_intervals.tsv` file they would not get in v0 without also passing `--confidence`. The extra file is informational and does not affect other outputs.
