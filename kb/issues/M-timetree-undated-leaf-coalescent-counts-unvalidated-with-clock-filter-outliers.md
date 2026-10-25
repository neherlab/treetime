# Undated-leaf lineage counts are not validated on runs with clock-filter outliers

Under a coalescent prior, every bad branch enters the expected lineage counts by the survival function of its posterior, and a leaf among them carries its leaf factor ([kb/decisions/timetree-undated-leaves-enter-lineage-counts-by-posterior.md](../decisions/timetree-undated-leaves-enter-lineage-counts-by-posterior.md)). A tip that the clock filter flags is a bad branch, like a tip without a date. The validation of the decision ran with `--clock-filter 0` (`validation/run.py`), so it covers undated tips but no outliers.

## Evidence

`dev/smoke`, the commit that counts bad branches by their posterior against `rust` before it, `timetree` cases with a coalescent prior:

- 13 cases on `ebola/20`, `flu/h3n2/20` and `zika/20`, where every leaf is dated and the clock filter flags no tip: byte-identical
- `timetree/ebola/362/skyline` (`--coalescent=skyline --seed=42`), where the clock filter flags 18 tips: the clock rate moves from $8.04 \cdot 10^{-4}$ to $7.93 \cdot 10^{-4}$ ($-1.4\%$), the skyline changes shape (segment boundaries move, $T_c$ rises up to 27-fold in the oldest segments and falls to 0.02 of before in the youngest; the first segment from 0.149 to 0.375 years), the coalescent log-likelihood falls from $-695.6$ to $-758.7$, and node dates move by up to 0.10 years

Whether counting a flagged tip by its posterior is closer to the truth than counting it as surviving at every time, as before, is not measured.

> [!IMPORTANT]
> **Investigation required.** Measure the effect against known truth: run the `opt-est` and `skyline-est` variants of the timetree suite on the messy validation data with the default clock filter (`--clock-filter 3`) instead of `--clock-filter 0`, before and after the decision, and compare the rate error, $\ln(T_c/T_c^\text{true})$ and the node date RMSE with the clean data. Record how many tips the filter flags per set.

## Validation

With the default clock filter on the messy data, the rate error, the $T_c$ error and the node date RMSE after the decision are no worse than before it, and close to those of the clean data.
