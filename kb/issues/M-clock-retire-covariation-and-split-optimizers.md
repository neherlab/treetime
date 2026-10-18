# Retire the covariation regression, its variance flags and the split optimizers

Phase 0c of [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md).
Decided: `--covariation` and `--tip-slack` are retired, not remapped.

## Scope

- `ClockVarianceParams` and the variance-carrying `ClockSet::propagate_averages`; the flags `--covariation`,
  `--tip-slack`, `--variance-factor`, `--variance-offset`, `--variance-offset-leaf`.
- `clock/find_best_root/` grid, Brent and golden-section split optimizers and `--branch-split-*`; the duplicate
  search in `reroot/` (the `optimize` min-dev search becomes the fixed-rate-0 case of the clock search).
- Outputs: the clock-model JSON statistics, the augur `rate_std` (from 0b), contracts, JSON schemas and the app's
  covariation switch, with the engineer.
- Tests: the clock and covariation snapshots of `dev/smoke.toml` are re-baselined with the reason recorded (the clock
  and reroot snapshots already change with 0a, which places the least-squares root exactly and keeps outliers out of
  the fit); the v0 runner golden masters move to the new API.

## kb

Resolves or obsoletes M-clock-covariation-variance-diverges-from-v0, M-clock-filter-regression-uses-covariation,
M-timetree-covariation-tip-slack-default-differs-from-v0, N-clock-regression-all-negative-rate,
N-reroot-clock-search-duplicates-generic-module, N-reroot-split-optimizer-default-diverges-from-v0,
N-reroot-leaf-variance-offset-diverges-from-v0. Needs a decision file for the replacement (v0 divergence).
