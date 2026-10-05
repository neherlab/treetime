# Clock filter residual computation differs from v0

> [!WARNING]
> **Needs review.** The dengue/100 measurement under Impact was not re-measured against the current code. It may include an interval-date effect: v1 now resolves interval dates to their mean, like v0 (see Differences).

v1's `clock_filter()` and v0's `residual_filter()` produce different outlier sets on the same data due to implementation differences in the IQD computation and outlier exclusion rules.

## Differences

1. IQD computation: v1 uses integer-rank indexing (`(3*n)/4` and `n/4`) for quartiles [packages/treetime/src/clock/clock_filter.rs#L60-L62](../../packages/treetime/src/clock/clock_filter.rs#L60-L62). v0 uses `np.percentile(residuals, 75)` which interpolates between adjacent values [packages/legacy/treetime/treetime/clock_filter_methods.py#L15](../../packages/legacy/treetime/treetime/clock_filter_methods.py#L15). For small datasets, this produces different IQD values.

2. Root-child exclusion: v0 skips children of the root from outlier flagging (`node.up.up is not None` check at [packages/legacy/treetime/treetime/clock_filter_methods.py#L18](../../packages/legacy/treetime/treetime/clock_filter_methods.py#L18)). v1 does not exclude any leaves based on tree position [packages/treetime/src/clock/clock_filter.rs#L64-L68](../../packages/treetime/src/clock/clock_filter.rs#L64-L68).

Date handling matches: v0 uses `np.mean(node.raw_date_constraint)` for interval dates, and v1's `likely_time()` returns the node time that `assign_dates` sets to the interval mean [packages/treetime/src/clock/clock_state.rs#L88-L90](../../packages/treetime/src/clock/clock_state.rs#L88-L90) [packages/treetime/src/clock/assign_dates.rs#L31](../../packages/treetime/src/clock/assign_dates.rs#L31).

## Impact

Different outlier sets lead to different final clock models. On dengue/100: v0 flags 8 outliers (R²=0.93), v1 flags 10 (R²=0.66). 7 of 8 v0 outliers overlap with v1's set.

## v0 Reference

[packages/legacy/treetime/treetime/clock_filter_methods.py#L5-L40](../../packages/legacy/treetime/treetime/clock_filter_methods.py#L5-L40) (`residual_filter`)

## v1 Location

[packages/treetime/src/clock/clock_filter.rs#L22](../../packages/treetime/src/clock/clock_filter.rs#L22) (`fn clock_filter`)
