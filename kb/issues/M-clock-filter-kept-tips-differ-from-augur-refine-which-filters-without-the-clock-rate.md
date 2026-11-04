# Kept tips differ from `augur refine`, which filters without the clock rate

## Problem

`augur refine` runs the clock filter of TreeTime v0 itself, right after it builds the `TreeTime` object and before the time-tree run, and passes no clock rate: `tt.clock_filter(reroot=reroot, n_iqd=clock_filter_iqd, plot=False)` in `refine()` of `augur/refine.py` (augur `58f32904`, line 54). With `--clock-rate`, the v0 filter therefore fits the slope of the root-to-tip regression, on the root-to-tip distances of the input branch lengths. v1 `timetree --clock-filter` filters with the given rate, after its branch-length optimization and its rerooting, as v0 does in its own `TreeTime.run()`, which passes `fixed_clock_rate` to the filter (`packages/legacy/treetime/treetime/treetime.py#L252-L258`).

When v1 replaces TreeTime v0 in `augur refine`, builds that fix the clock rate keep different tips, and the dates of the deep nodes then compare trees with different tips.

## Evidence

The comparison of TreeTime v0 and v1 in the Nextstrain workflows (pipelandia, round 3; v1 at `0e522d97e`, v0 through `augur refine` 34.3.0 with TreeTime 0.12.1, same commands, inputs, and image) finds 15 trees with tips that only one engine keeps. 12 of the 13 published trees among them fix the clock rate.

- **Threshold**: for every tip that only one engine keeps, the residuals of the two engines lie on opposite sides of the threshold
- **Refit with the fixed rate**: a regression of the v0 root-to-tip distances with the slope fixed at the clock rate of the build flags exactly the outliers of v1 in 7 trees, and those of v1 except 1 or 2 tips near the threshold in norovirus GII.17/genome, rubella genome, seasonal-cov hku1, and seasonal-cov oc43
- **seasonal-cov hku1**: the v0 filter fits a slope of -8.9e-05 substitutions per site per year, against the fixed rate of 0.00012. v0 sees no clock signal and flags no tip, while v1 flags 18
- **Other causes**: in norovirus all/RdRp, v1 roots the tree next to the 3 tips that v0 flags, so their v1 residuals are 0.05 to 0.11 interquartile distances; seasonal-cov nl63 (`--root best`) and mumps north-america (estimated rate) each differ by one tip near the threshold
- **v1 flagging rule**: v1 flags outliers with the quartiles and the root rule of v0 (`0990e2b04`, `35a713a94`); this change left the kept tips of the published builds unchanged

> [!IMPORTANT]
> **Decision required.** The augur integration of v1 must choose one filter. Options:
>
> - **Reproduce the augur filter**: v1 offers the filter that `augur refine` runs (the regression slope also when the build fixes the rate, on the input branch lengths, before branch-length optimization and rerooting), so `augur refine` with v1 keeps the tips that it keeps with v0
> - **Replace the augur filter**: `augur refine` with v1 uses the v1 filter with the fixed rate, and the integration records the changed kept tips as an intended difference
> - **Change augur**: `augur refine` passes the clock rate of the build to the filter, so the filter of v0 through augur and the filter of v1 use the same rate
>
> Evidence for the choice: on seasonal-cov hku1 the regression slope is negative, so the augur filter has no clock signal to filter by, while the fixed rate of the build gives one.

## Reproduction

The residuals of both filters come from the `augur refine` command of each build: v0 with a wrapper around `residual_filter()` of TreeTime v0 that records the clock model, the interquartile distance, and the residual of every dated tip before the filter marks outliers; v1 through the adapter of the pipelandia comparison with the clock-filter threshold lowered to 1e-12 and without pruning, so that its outlier table lists every dated tip that is not a child of the root.
