# v1 clock regression produces all-negative rates where v0 finds positive

On some datasets (e.g. dengue/100 with outliers), v1's clock regression estimates negative clock rate at all root positions during the pre-filter step, while v0 finds at least some positive-rate positions using the same `force_positive=True` constraint.

`treetime clock --covariation` on `data/flu/h3n2/20` (with `--sequence-length=1400`) fails the same way: "the pre-filter step removed outliers but the clock rate remains negative at all root positions". Without `--covariation` the command succeeds. v0 on the same input with `--covariation` estimates a rate of 2.765e-03 ± 4.03e-04 with R² 0.98 and a root date of 1996.72. This is the `clock/flu/h3n2/20/covariation` case of `dev/smoke`, declared an expected failure in `dev/smoke.toml`; the `rust` branch fails it too.

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
