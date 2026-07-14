# Wire overdispersion and --tip-slack through to clock regression

v1 hardcodes `overdispersion = 2.0` for the covariation clock path. v0 uses `OVER_DISPERSION = 10` (stored as `tip_slack`), with CLI override via `--tip-slack`.

- v1: `overdispersion = 2.0` at
  [`packages/treetime/src/commands/clock/run.rs#L64`](../../packages/treetime/src/commands/clock/run.rs#L64),
  fed into `ClockParams.variance_factor` as `overdispersion / seq_len`
- v0: `self.tip_slack = ttconf.OVER_DISPERSION` at
  [`packages/legacy/treetime/treetime/clock_tree.py#L98`](../../packages/legacy/treetime/treetime/clock_tree.py#L98),
  constant `OVER_DISPERSION = 10` at [`packages/legacy/treetime/treetime/config.py#L8`](../../packages/legacy/treetime/treetime/config.py#L8), used in branch variance at [`packages/legacy/treetime/treetime/clock_tree.py#L278-L285`](../../packages/legacy/treetime/treetime/clock_tree.py#L278-L285)

## Background

Molecular evolution shows overdispersion relative to a Poisson process: observed branch-length variance exceeds the Poisson expectation (mean = variance). The overdispersion parameter inflates variance estimates for terminal branches in the weighted least squares clock regression, giving outlier tips less influence on the fitted rate.

v0 adds `tip_slack^2 * one_mutation^2` to terminal branch variance, where `one_mutation = 1 / seq_len`. With `tip_slack = 10` and `seq_len = 1000`, the extra variance is `(10/1000)^2 = 1e-4`.

v1 uses `variance_factor = overdispersion / seq_len`. With `overdispersion = 2` and `seq_len = 1000`, the factor is `2e-3`, applied differently in the variance formula.

## Impact

The smaller overdispersion value in v1 reduces slack for terminal branches, making the clock regression more sensitive to outlier tips when `--covariation` is enabled. The v1 TODO comment flags this as provisional.

The `--tip-slack` CLI argument is not wired through in v1, so users cannot override the hardcoded value.

## Fix

Wire `--tip-slack` through to `ClockParams.variance_factor` and match v0's default `OVER_DISPERSION = 10`.

## Related issues

- Source: [kb/issues/M-clock-covariation-overdispersion.md](../issues/M-clock-covariation-overdispersion.md) -- delete after full resolution
