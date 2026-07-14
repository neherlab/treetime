# Clock filter residual computation differs from v0

v1's `clock_filter_inplace()` and v0's `residual_filter()` produce different outlier sets on the same data due to implementation differences in the IQD computation and outlier exclusion rules.

## Differences

1. IQD computation: v1 uses integer-rank indexing (`(3*n)/4` and `n/4`) for quartiles. v0 uses `np.percentile(residuals, 75)` which interpolates between adjacent values. For small datasets, this produces different IQD values.

2. Root-child exclusion: v0 skips children of the root from outlier flagging (`node.up.up is not None` check at `clock_filter_methods.py:18`). v1 does not exclude any leaves based on tree position.

3. Date handling: v0 uses `np.mean(node.raw_date_constraint)` for interval dates. v1 uses `likely_time()` which may resolve intervals differently.

## Impact

Different outlier sets lead to different final clock models. On dengue/100: v0 flags 8 outliers (R²=0.93), v1 flags 10 (R²=0.66). 7 of 8 v0 outliers overlap with v1's set.

## v0 Reference

`packages/legacy/treetime/treetime/clock_filter_methods.py:1-36` (`residual_filter`)

## v1 Location

`packages/treetime/src/clock/clock_filter.rs:21-105` (`clock_filter_inplace`)

## Related tickets

- [kb/tickets/clock-align-filter-residual-with-v0.md](../tickets/clock-align-filter-residual-with-v0.md)
