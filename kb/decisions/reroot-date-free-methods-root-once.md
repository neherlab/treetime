# Date-free rooting methods root once, on the starting branch lengths

The date-free rooting methods `min-dev`, `--outgroup` and `midpoint` choose the root from the tree alone, in every command. Each of them roots once, on the starting branch lengths, and counts every tip. The root then stays for the rest of the run, and every clock fit of the run still happens. The date-based methods `least-squares` and `oldest` keep their schedule.

**Type**: Scientific design choice (differs from v0). Approved by the user on 2026-10-24.

**v1 location**: `fn reroot_date_free()` in `packages/treetime/src/reroot/date_free.rs`, called by `fn estimate_clock_model_with_prefilter()` (`packages/treetime/src/clock/pipeline.rs`), `fn reroot_and_fit_to_dates()` (`packages/treetime/src/timetree/optimization/reroot.rs`) and `fn reroot_optimize()` (`packages/treetime/src/optimize/pipeline.rs`). `fn RerootSpec::uses_dates()` (`packages/treetime/src/reroot/params.rs`) separates the two groups.

## Decision

- **Starting branch lengths**: the input lengths in `clock` and `timetree`. In `optimize`, the lengths after the initial guess, which follows `--initial-branch-lengths`
- **One reroot**: `timetree` reroots before its first clock fit and skips the two reroots of the pre-loop, while both clock refits of the pre-loop stay. `clock` reroots before its first clock fit, and the fit after the outlier filter keeps the root. `optimize` reroots right after the initial guess, and its main loop optimizes the branch lengths on the rooted tree
- **Every tip counts**: the root search gives every tip the same weight, whether it has a date or not, and leaves no clock-filter outlier out. A dated stem root (a root with one child) stays in `clock` and `timetree`, and becomes a tip; `optimize` has no dates and removes it
- **Valid lengths**: `min-dev` and `midpoint` place the root from branch lengths and stop with an input error when an edge has a missing, negative or non-finite length. `--outgroup` needs no lengths

## Reasons

- A user picks a date-free method so that the dates do not decide the root, for example when the temporal signal is too weak to root by regression. The later reroots of `timetree` and `clock` leave out the outliers that the clock filter finds from the dates, so a reroot there would let the dates move the root
- The later branch lengths are not better estimates for the root search: the timetree pre-step gives 75% of the input lengths and 25% of one optimization pass, and converged lengths exist only after the main loop
- Under a reversible substitution model the likelihood does not depend on the root position (the pulley principle), so all later steps can run on the final root
- The three commands give the same root on the same tree with the same starting lengths, and `midpoint` gives the root that augur and Biopython give

## v0 divergence

TreeTime v0 reroots at several points of a run, with every method: in `TreeTime._run()` and in `TreeTime.clock_filter()` [packages/legacy/treetime/treetime/treetime.py#L260](../../packages/legacy/treetime/treetime/treetime.py#L260), [#L279](../../packages/legacy/treetime/treetime/treetime.py#L279), [#L487](../../packages/legacy/treetime/treetime/treetime.py#L487) and [#L513](../../packages/legacy/treetime/treetime/treetime.py#L513). Its `min_dev` score counts only the dated tips that are not outliers. v0 has no `midpoint` method, and `augur refine --root mid_point` roots with Biopython before TreeTime runs.

## Outputs this changes

- `min-dev` and `--outgroup` in `timetree`, and in `clock` with the outlier filter
- `min-dev` in `clock` and `timetree` when some tips have no date
- every `optimize` run that reroots, because the damped branch-length pass before the reroot is gone
- `min-dev` on a tree with a missing, negative or non-finite branch length now stops with an input error; before, the root search counted a missing length as zero and used a negative length as it was

## Related

- [kb/decisions/clock-min-dev-fixed-zero-rate.md](clock-min-dev-fixed-zero-rate.md): the `min-dev` score
- [kb/decisions/reroot-named-node-at-unrooted-branch-middle.md](reroot-named-node-at-unrooted-branch-middle.md): the root that `--outgroup` places
- [kb/issues/M-timetree-reroot-schedule-diverges-from-v0.md](../issues/M-timetree-reroot-schedule-diverges-from-v0.md): the open decision on the schedule of `least-squares` and `oldest`
