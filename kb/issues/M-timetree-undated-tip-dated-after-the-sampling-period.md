# An undated tip on a long branch is dated after the sampling period

A tip without any date (metadata `XXXX-XX-XX`) has no date constraint. The forward pass places it by the message from its parent alone: its time posterior is the parent's posterior convolved with the branch likelihood of its divergence (`fn refine_distribution_from_parent()` in [packages/treetime/src/timetree/inference/forward_pass.rs](../../packages/treetime/src/timetree/inference/forward_pass.rs)). Nothing bounds this posterior from above, so the tip lands at its parent's date plus the time that the clock needs for the divergence of its branch:

$$t_{\text{tip}} \approx t_{\text{parent}} + \frac{d}{\mu}$$

where:

- $t_{\text{parent}}$ -- date of the parent node
- $d$ -- divergence of the tip's branch in substitutions per site (`mutation_length`)
- $\mu$ -- clock rate in substitutions per site per year

When $d$ is large, this date lies after the latest sample, or after the present. The run gives no warning: `fn dates_outside_input_dates()` [packages/treetime/src/timetree/dates_outside_input.rs#L11](../../packages/treetime/src/timetree/dates_outside_input.rs#L11) checks tips against their input dates only and skips tips that have none.

## Evidence

lassa l production case of pipelandia (752 tips, 18 with `XXXX-XX-XX`, 228 with a year only, latest full date 2024-09-21), `--clock-rate 0.00078 --keep-root --confidence`, v1 at `rust` `34055d3bd` against augur 34.3.0 with TreeTime 0.12.1 (`augur refine --keep-root --date-inference marginal --date-confidence`) on the same tree and alignment:

| tip | coalescent | v1 date [90%] | v0 date [90%] | parent date v1 / v0 | $d$ |
| --- | --- | --- | --- | --- | --- |
| FV537337 | `opt` | 2254.41 [2241.40, 2267.24] | 2254.19 [2241.43, 2267.15] | 1979.54 / 1979.58 | 0.2142 |
| FV537337 | none | 2254.40 [2241.39, 2267.24] | 2254.19 [2241.42, 2267.14] | 1979.54 / 1979.58 | 0.2142 |
| U73035 | none | 2020.23 [2015.78, 2024.99] | 2020.09 [2015.77, 2024.98] | 1979.12 / 1979.16 | 0.0319 |

- $d/\mu = 0.2142 / 0.00078 = 274.6$ years, and v1 reports a branch of 274.9 years, so the date of `FV537337` is the model's answer. The width of its interval matches the Poisson spread of about 1559 substitutions on 7279 sites, $\sqrt{dL}/(L\mu) \approx 7$ years of standard deviation, plus the parent's uncertainty
- The other 17 undated tips lie within the sampling period, and v1 and v0 agree within 1 year on each of them except `AY179171` and `FR832710`, which sit on zero-length branches below a deep node whose date differs between the engines by about 35 years (v1 1896.0, v0 1930.6 without the coalescent prior). On `OR147788`, also on a zero-length branch, the dates agree, but v0 reports the interval [2009.48, 2041.80] against [2017.84, 2018.47] in v1. On the year-only tips, the year range bounds the date, and no tip leaves its year in either engine
- The coalescent prior does not move these dates: an undated leaf enters the lineage counts by its posterior ([kb/decisions/timetree-undated-leaves-enter-lineage-counts-by-posterior.md](../decisions/timetree-undated-leaves-enter-lineage-counts-by-posterior.md)), and past the latest sample its leaf factor is nearly flat

## v0 behavior

v1 follows v0. v0 propagates the parent message to every node without a delta constraint, unconstrained terminals included ([packages/legacy/treetime/treetime/clock_tree.py#L790-L860](../../packages/legacy/treetime/treetime/clock_tree.py#L790-L860)), and augur writes the resulting `numdate` with its confidence interval for undated tips.

## Impact

- `timetree.node-data.json`, `timetree.auspice.json`, and the Newick and Nexus outputs carry a tip date centuries after the samples, with `date_inferred: true`, and a branch of the same length; the time axis of an Auspice view stretches to that date. `timetree.dates.tsv` writes `--` for undated tips
- The coalescent lineage grid reaches past the latest sample for the same cause ([M-coalescent-lineage-grid-extends-past-the-latest-sample-with-undated-tips.md](M-coalescent-lineage-grid-extends-past-the-latest-sample-with-undated-tips.md))
- `FV537337` is a patent record of GenBank, a kind of entry that often carries no collection date

## Fix direction

> [!IMPORTANT]
> **Decision required.** The date is what v0 and the model give, so each option except (a) changes scientific outputs or adds a divergence from v0. Options:
> - (a) keep the v0 behavior, and extend the run warning `dates-outside-input-dates` to undated tips dated after the latest input date, so the user sees the case
> - (b) bound the posterior of an undated tip by a sampling window, such as the range of the observed dates, in the forward pass only, so $\log Z$ stays free of a sampling-date prior; this is option (b) of [M-coalescent-lineage-grid-extends-past-the-latest-sample-with-undated-tips.md](M-coalescent-lineage-grid-extends-past-the-latest-sample-with-undated-tips.md), and one choice settles both issues
> - (c) report no date for an undated tip in the tree outputs, as `timetree.dates.tsv` already does, and keep its posterior for the coalescent and the parent's date
>
> Evidence: the table above. Option (b) would place `FV537337` at the end of the window and leave the other 17 undated tips of lassa l nearly in place, because their dates lie inside it (not measured).
