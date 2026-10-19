# Clock-filter outliers lose their date constraint, so outlier clades can be dated after the samples

The clock filter marks tips whose root-to-tip divergence does not fit the clock as bad branches
(`fn bad_leaves()` [packages/treetime/src/timetree/inference/bad_branches.rs#L7](../../packages/treetime/src/timetree/inference/bad_branches.rs#L7)).
A bad branch sends no message to its parent in the backward pass, so the time inference dates a flagged tip and the
ancestors above flagged tips from the clock alone. When a whole clade is flagged, the clock can place the tip dates
outside their input date ranges and the internal nodes of the clade after the latest sample. Such an internal node is
dated after every one of its descendants.

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

- `data/mpox/clade-ii/1000` with `--clock-rate=5.7e-05` (default flags otherwise): the clock filter flags 112 tips.
  42 tips lie outside their input date ranges, by up to 9.52 years; 91 internal nodes lie after the latest input
  date 2024.62, by up to 5.63 years (the latest at 2030.25, above two tips sampled in 2023). `rust` flags 8 tips on
  the same run; `rust` with `--keep-root` flags 98 tips and gives the same dates
- `data/ebola/362` with default flags: the clock filter flags 21 tips. 10 internal nodes of one clade of `EMLH`
  samples (sampled 2015.4 to 2015.6) lie after the latest input date 2016.25, by up to 0.96 years. The run also
  reports 14 branches that run backwards in time
- `data/flu/h3n2/200` with default flags: no date outside the input dates

## v0 behavior

v0 also ignores the date constraint of a bad branch: the backward pass gives a bad branch no information
([packages/legacy/treetime/treetime/clock_tree.py#L685-L687](../../packages/legacy/treetime/treetime/clock_tree.py#L685-L687)),
the forward pass keeps a date constraint only on nodes that are not bad branches
([packages/legacy/treetime/treetime/clock_tree.py#L792](../../packages/legacy/treetime/treetime/clock_tree.py#L792)),
and the setup logs that the date constraint of a bad branch will be ignored
([packages/legacy/treetime/treetime/clock_tree.py#L389-L396](../../packages/legacy/treetime/treetime/clock_tree.py#L389-L396)).

v1 differs for flagged tips with an exact date: the forward pass leaves a node with an exact date unrefined
(`fn refine_distribution_from_parent()` [packages/treetime/src/timetree/inference/forward_pass.rs#L137](../../packages/treetime/src/timetree/inference/forward_pass.rs#L137)),
so such a tip keeps its input date while its ancestors are dated from the clock. This gives the branches that run
backwards in time on `data/ebola/362`. Flagged tips with a date range are dated from the clock, as in v0.

> [!IMPORTANT]
> **Decision required.** The model fix is open. Options:
>
> - Keep the v0 behavior and the run warning: flagged tips and their ancestors stay dated by the clock alone, and the
>   warning names them
> - Keep the date range of flagged tips in the time inference: a flagged tip leaves the clock fit but still bounds
>   its own date and the dates of its ancestors
> - A local clock filter: flag a tip by its divergence relative to its neighbors in the tree instead of the global
>   root-to-tip regression, so that a clade with its own rate is not flagged as a whole
>
> Evidence to weigh: the counts above, and the date and rate errors of each option on the validation suite.

## Validation

On `data/mpox/clade-ii/1000` with `--clock-rate=5.7e-05` and on `data/ebola/362` with default flags, no tip lies
outside its input date range and no internal node lies after the latest sample, or the warning names each one that
the chosen option accepts.
