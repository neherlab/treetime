# The coalescent lineage grid extends past the latest sample with undated tips

`fn LineageExpectation::span()` [packages/treetime/src/coalescent/mean_field.rs#L467-L477](../../packages/treetime/src/coalescent/mean_field.rs#L467-L477) ends at the last point of the expectation grid, which reaches the latest time at which any node's survival function is below one. An undated leaf enters the lineage counts by its posterior ([kb/decisions/timetree-undated-leaves-enter-lineage-counts-by-posterior.md](../decisions/timetree-undated-leaves-enter-lineage-counts-by-posterior.md)). Its posterior has no upper bound from a sampling window: past the latest sample $\bar k_{-v} \approx 0$, so its leaf factor is nearly flat there, and the right tail follows the branch likelihood.

## Evidence

- **Validation**: `validation/data_n100_messy` (3% of tips undated), `keep-opt-fixed`, $L=10^3$, $\mu=3\cdot10^{-4}$, `const`: the one segment of the coalescent output ends between 2057 and 2088 in the 10 replicates, against a latest sample date of 2025.08. $T_c$ is 17.0 against 16.5 on the clean data, so the fit itself is barely affected
- **Evaluation**: `validation/evaluate.py` scores a segment at its midpoint and drops segments whose midpoint lies after the present, so `tc_log_bias` of that cell is `NaN`

## Impact

- `fn optimize_skyline()` [packages/treetime/src/coalescent/skyline.rs#L40-L41](../../packages/treetime/src/coalescent/skyline.rs#L40-L41) places equal-width segments over the span, so segments after the latest sample hold no lineages and the segments that hold data are wider
- The coalescent outputs (`--output-coalescent`) report a segment end decades after the samples

> [!IMPORTANT]
> **Investigation required.** The effect on skyline $T_c$ by depth and on its band coverage on the messy validation data at low signal, against the clean data.

## Fix direction

> [!IMPORTANT]
> **Decision required.** Options: (a) end the span at the latest dated leaf, or at the latest time where $\bar k \ge 1$, and let the survival tails beyond it fall in the last segment; (b) bound the posterior of an undated leaf by a sampling window, such as the range of the observed dates, in the forward pass only, so $\log Z$ stays free of a sampling-date prior. Option (a) changes the segments only; option (b) also changes the dates of undated leaves.
