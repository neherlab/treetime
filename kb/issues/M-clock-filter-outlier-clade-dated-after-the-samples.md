# Clock-filter outliers lose their date constraint, so outlier clades can be dated after the samples

The clock filter marks tips whose root-to-tip divergence does not fit the clock as bad branches
(`fn bad_leaves()` [packages/treetime/src/timetree/inference/bad_branches.rs#L11](../../packages/treetime/src/timetree/inference/bad_branches.rs#L11)).
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

> [!IMPORTANT]
> **Decision required.** The default run follows v0, which ignores the date constraint of a bad branch. Any change of
> the default dates differs from v0 and needs user approval. Options: (a) keep the v0 behavior and keep only the
> `dates-outside-input-dates` warning; (b) make `--prune-filtered-tips` the default for `timetree` and `clock`, as augur
> refine does; (c) keep the flagged tips but apply their input date as a soft constraint (a date range widened by the
> clock-filter residual) in both time passes, so no flagged node is dated from the clock alone.

## Proposed implementation plan

The current code matches the description: `fn bad_leaves()` and `fn applied_date_constraint()` in
`packages/treetime/src/timetree/inference/bad_branches.rs`, used by `forward_pass.rs` and `backward_pass.rs` of the same
directory; `fn dates_outside_input_dates()` in `packages/treetime/src/timetree/dates_outside_input.rs`; pruning in
`fn prune_tips()` of `packages/treetime/src/timetree/pipeline.rs` and `packages/treetime/src/clock/prune_tips.rs`,
behind `prune_filtered_tips` of `packages/app-commands/src/commands/shared/clock_filter.rs` (default off).

1. Get the user's choice of the decision above, and record an approved divergence from v0 in `kb/decisions/`
2. Option (a): no code change; delete this issue after a note in the `--help` text of `--clock-filter` that flagged tips
   keep only clock dates unless `--prune-filtered-tips` is given
3. Option (b): change the default of `prune_filtered_tips` (with a negation flag, per the clap rules) in
   `clock_filter.rs`, update `packages/app-cli/src/cli/v0_args.rs` only if a v0 spelling changes meaning, run
   `just gen` for the schemas, and update `test_dependent_settings.rs` and `test_pipeline_prune.rs`
4. Option (c): extend `fn applied_date_constraint()` to return a widened constraint for a bad branch instead of none,
   with the width taken from the clock-filter residual; both passes read it unchanged
5. Tests: a unit test in `packages/treetime/src/timetree/inference/__tests__/` on a small hand-made tree with one
   flagged clade, asserting that no internal node is dated after the latest sample (oracle: the definition of the
   input date range); for option (c) also an exact-date flagged tip whose posterior lies within the widened range
   (oracle: manual calculation of the constrained posterior on the small tree)
6. Validation: the runs of the Validation section above with `just r treetime ...`, checking that the
   `dates-outside-input-dates` warning is absent; `just check`, `just test-rs`, and `just generated-check` for option (b)
