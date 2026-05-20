# v1 clock regression produces all-negative rates where v0 finds positive

On some datasets (e.g. dengue/100 with outliers), v1's clock regression estimates negative clock rate at all root positions during the pre-filter step, while v0 finds at least some positive-rate positions using the same `force_positive=True` constraint.

## Impact

The root cause is not fully understood. v1 compensates by using `force_positive_rate: false` for the pre-filter step (documented in `decisions/clock-prefilter-relaxed-positive-rate.md`). This works but means v1's pre-filter uses a different root than v0 would choose, contributing to different outlier sets.

## Possible causes

- Different variance model defaults between v0 and v1
- Differences in how branch lengths or dates are processed before regression
- Numerical differences in the regression accumulation order

## v0 Reference

v0's pre-filter reroot at `packages/legacy/treetime/treetime/treetime.py:486-489` succeeds with `force_positive=True` on dengue/100.

## v1 Location

`packages/treetime/src/clock/find_best_root/find_best_root.rs`: all 198 nodes rejected for negative rate on dengue/100 pre-filter.

## Related issues

- Source: [N-clock-regression-all-negative-rate.md](../issues/N-clock-regression-all-negative-rate.md) -- delete after full resolution
- [clock-prefilter-relaxed-positive-rate.md](../decisions/clock-prefilter-relaxed-positive-rate.md) -- v1 workaround for this issue
