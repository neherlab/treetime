# The coalescent lineage grid extends past the latest sample with undated tips

`fn LineageExpectation::span()` [packages/treetime/src/coalescent/mean_field.rs#L457-L467](../../packages/treetime/src/coalescent/mean_field.rs#L457-L467) ends at the last point of the expectation grid, which reaches the latest time at which any node's survival function is below one. An undated leaf enters the lineage counts by its posterior ([kb/decisions/timetree-undated-leaves-enter-lineage-counts-by-posterior.md](../decisions/timetree-undated-leaves-enter-lineage-counts-by-posterior.md)). Its posterior has no upper bound from a sampling window: past the latest sample $\bar k_{-v} \approx 0$, so its leaf factor is nearly flat there, and the right tail follows the branch likelihood.

Clock-filter outliers extend the span the same way. A tip that the clock filter flags loses its date constraint, and the time inference dates it from the clock alone ([M-clock-filter-outlier-clade-dated-after-the-samples.md](M-clock-filter-outlier-clade-dated-after-the-samples.md)), so its posterior has no upper bound from the sampling window either.

## Evidence

- **Validation**: `validation/data_n100_messy` (3% of tips undated), `keep-opt-fixed`, $L=10^3$, $\mu=3\cdot10^{-4}$, `const`: the one segment of the coalescent output ends between 2057 and 2088 in the 10 replicates, against a latest sample date of 2025.08. $T_c$ is 17.0 against 16.5 on the clean data, so the fit itself is barely affected
- **Clock-filter outliers**: `data/ebola/362`, `--coalescent skyline --n-skyline 50 --clock-filter 4` (a web app run): the 344 tips that the clock model keeps lie between 2014.21 and 2015.81, the 18 flagged tips are dated between 2015.08 and 2016.35, and the skyline segments end at 2018.66, so 30 of the 50 segments start after the latest kept sample. The app plots the skyline up to the latest kept sample (`latest_sample` of `TimetreeEstimates` in [packages/app-backend/src/results/timetree.rs](../../packages/app-backend/src/results/timetree.rs)); the coalescent output still reports every segment
- **Evaluation**: `validation/evaluate.py` scores a segment at its midpoint and drops segments whose midpoint lies after the present, so `tc_log_bias` of that cell is `NaN`

## Impact

- `fn optimize_skyline()` [packages/treetime/src/coalescent/skyline.rs#L40-L41](../../packages/treetime/src/coalescent/skyline.rs#L40-L41) places equal-width segments over the span, so segments after the latest sample hold no lineages and the segments that hold data are wider
- The coalescent outputs (`--output-coalescent`) report a segment end decades after the samples

> [!IMPORTANT]
> **Investigation required.** The effect on skyline $T_c$ by depth and on its band coverage on the messy validation data at low signal, against the clean data.

## Fix direction

> [!IMPORTANT]
> **Decision required.** Options: (a) end the span at the latest dated leaf, or at the latest time where $\bar k \ge 1$, and let the survival tails beyond it fall in the last segment; (b) bound the posterior of an undated leaf by a sampling window, such as the range of the observed dates, in the forward pass only, so $\log Z$ stays free of a sampling-date prior. Option (a) changes the segments only; option (b) also changes the dates of undated leaves, and leaves the span of clock-filter outliers unbounded unless it bounds their posterior too. For option (a) the latest dated leaf is the latest leaf whose date the clock model keeps.

## Proposed implementation plan

The locations above are current: `fn LineageExpectation::span()` takes `t_max` from the last grid point (`mean_field.rs` L465), and `fn optimize_skyline()` builds equal-width boundaries over that span (`skyline.rs` L40-L41). The plan depends on the open decision.

- **Option (a)**:
  1. Pass the latest kept dated leaf (the latest leaf whose date the clock model keeps) into `LineageExpectation`, and end `span()` at it, or at the latest grid time with $\bar k \ge 1$
  2. In `fn optimize_skyline()`, let the last segment take the survival mass beyond `t_max`, so the pair integral and the merger counts of the tails stay in the fit
  3. Tests in `packages/treetime/src/coalescent/__tests__/`: a small tree with one undated leaf whose posterior reaches past the latest sample; assert that the span ends at the latest dated leaf and that the sum of the segment mergers equals the total of the unbounded grid (oracle: conservation of the expected merger count)
- **Option (b)**:
  1. Record the approved choice in `kb/decisions/` next to `timetree-undated-leaves-enter-lineage-counts-by-posterior.md`
  2. In the forward pass, multiply the posterior of an undated leaf (and of a clock-filter outlier, if approved) by an indicator of the range of the observed dates, and keep the backward messages and $\log Z$ unchanged
  3. Tests: an undated leaf on a small tree whose posterior support ends at the latest observed date (oracle: the definition of the sampling window), and an unchanged $\log Z$ against the current code (oracle: $\log Z$ comes from the backward pass, which the option does not change)
- **Validation for both**: `just check`, `just test-rs`; rerun `validation/data_n100_messy` (`keep-opt-fixed`, $L=10^3$, $\mu=3\cdot10^{-4}$, `const`) and `data/ebola/362` with `--coalescent skyline --n-skyline 50 --clock-filter 4`, and check that no segment ends after the latest kept sample and that `tc_log_bias` is finite. This also answers the Investigation required block above
