# Clock filter regression uses covariation

Before it detects outliers, the clock pre-filter fits a root-to-tip regression. v1 passes the command's variance parameters to that fit, so with `--covariation` the outlier filter runs on a covariation-aware regression (`fn estimate_clock_model_with_prefilter`, [packages/treetime/src/clock/pipeline.rs#L117-L160](../../packages/treetime/src/clock/pipeline.rs#L116-L160), argument `options`).

v0 always fits the filter regression without covariation: `TreeTime.clock_filter` calls `self.reroot(..., covariation=False, ...)` or `self.get_clock_model(covariation=False, ...)` ([packages/legacy/treetime/treetime/treetime.py#L457-L491](../../packages/legacy/treetime/treetime/treetime.py#L457-L491)).

## Impact

With `--covariation`, v1 can flag a different set of outlier tips than v0. The residuals that the filter thresholds come from a regression with different weights and possibly a different root. Without `--covariation`, both versions use the same unweighted-tip regression.

## Related

- [M-clock-covariation-variance-diverges-from-v0.md](M-clock-covariation-variance-diverges-from-v0.md)
- [M-clock-filter-residual-parity.md](M-clock-filter-residual-parity.md): other filter residual differences
