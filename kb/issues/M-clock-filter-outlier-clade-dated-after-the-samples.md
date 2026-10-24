# Clock-filter outliers lose their date constraint, so outlier clades can be dated after the samples

The clock filter marks tips whose root-to-tip divergence does not fit the clock as bad branches
(`fn bad_leaves()` [packages/treetime/src/timetree/inference/bad_branches.rs#L7](../../packages/treetime/src/timetree/inference/bad_branches.rs#L7)).
A bad branch sends no message to its parent in the backward pass, and neither time pass applies its input date, so
the time inference dates a flagged tip and the ancestors above flagged tips from the clock alone. The clock places
the flagged tips outside their input date ranges, and when a whole clade is flagged, the internal nodes of the clade
after the latest sample. Such an internal node is dated after every one of its descendants' input dates.

The run warning `dates-outside-input-dates` reports both cases after the final time inference
(`fn dates_outside_input_dates()` [packages/treetime/src/timetree/dates_outside_input.rs#L11](../../packages/treetime/src/timetree/dates_outside_input.rs#L11)).
For a tip with input date range $[a, b]$ and final date $t$, the warning reports the excess

$$e_{\text{tip}} = \max(a - t,\ t - b) > 0$$

and for an internal node with final date $t$ it reports

$$e_{\text{internal}} = t - b_{\max} > 0$$

where:

- $a$, $b$ -- lower and upper bound of the input date of the tip ($a = b$ for an exact date)
- $t$ -- final date of the node
- $b_{\max}$ -- the largest upper bound of all input dates of the run

## Evidence

Default flags unless noted, `--seed 42`:

- `data/sc2/4500` with the example config `data/sc2/4500/timetree.yaml` (clock filter 3, skyline coalescent): the clock filter flags 500 tips. All 500 lie outside their input dates, by up to 5.06 years, and 105 internal nodes lie after the latest input date 2026.36, by up to 1.64 years. The rate step of the coalescent prior lowers the clock rate to 0.72 times the regression rate, which moves these nodes later ([M-timetree-rate-step-sets-sc2-clock-rate-below-the-regression.md](M-timetree-rate-step-sets-sc2-clock-rate-below-the-regression.md))
- `data/sc2/4500`: the clock filter flags 410 tips; all 410 lie outside their input dates, by up to 5.06 years. No internal node lies after the latest input date
- `data/ebola/362`: the clock filter flags 18 tips; all 18 lie outside their input dates, by up to 1.06 years. No internal node lies after the latest input date
- `data/mpox/clade-ii/1000` with `--clock-rate=5.7e-05` (measured before flagged tips with an exact date were dated from the tree): the clock filter flags 112 tips; 91 internal nodes lie after the latest input date 2024.62, by up to 5.63 years (the latest at 2030.25, above two tips sampled in 2023)
- `data/flu/h3n2/200`: no date outside the input dates

## v0 behavior

v1 follows v0. v0 ignores the date constraint of a bad branch: the backward pass gives a bad branch no information
([packages/legacy/treetime/treetime/clock_tree.py#L685-L687](../../packages/legacy/treetime/treetime/clock_tree.py#L685-L687)),
the forward pass keeps a date constraint only on nodes that are not bad branches
([packages/legacy/treetime/treetime/clock_tree.py#L792](../../packages/legacy/treetime/treetime/clock_tree.py#L792)),
and the setup logs that the date constraint of a bad branch will be ignored
([packages/legacy/treetime/treetime/clock_tree.py#L389-L396](../../packages/legacy/treetime/treetime/clock_tree.py#L389-L396)).
In v1 both time passes read the date of a node through `fn applied_date_constraint()`
([packages/treetime/src/timetree/inference/bad_branches.rs](../../packages/treetime/src/timetree/inference/bad_branches.rs)),
which gives a bad branch no date, exact or range. A flagged tip is therefore dated by the message from its parent
alone and is never dated before its parent.

augur removes the flagged tips from the tree before the time inference
(`augur/refine.py` lines 52-62 of `nextstrain/augur` `58f32904`), so its outputs hold neither the flagged tips nor
the dates the clock gives them.

## Fix direction

A run of `timetree` or `clock` with `--prune-filtered-tips` removes the flagged tips and their clades from the tree
before the time inference, as augur refine does, so its outputs hold no date that comes from the clock alone.

The default run keeps the flagged tips, so it still dates outlier clades from the clock alone. This is the open part of
the issue.

## Validation

On `data/sc2/4500` with `data/sc2/4500/timetree.yaml`, on `data/ebola/362` with default flags, and on
`data/mpox/clade-ii/1000` with `--clock-rate=5.7e-05`, no tip lies outside its input date range and no internal node
lies after the latest sample, both with `--prune-filtered-tips` and with the default flags.
