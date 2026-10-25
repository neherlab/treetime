# Clock-filter outliers cannot be pruned before the time inference

## Problem

`augur refine --clock-filter-iqd` removes every tip that the TreeTime clock filter flags from the tree before the time inference, except the tips that `--keep-ids` lists in a file, one name per line ([augur/refine.py#L49-L62](https://github.com/nextstrain/augur/blob/58f329042e7fc631ca275cf2df44c8720f80949f/augur/refine.py#L49-L62), [#L213](https://github.com/nextstrain/augur/blob/58f329042e7fc631ca275cf2df44c8720f80949f/augur/refine.py#L213)). TreeTime itself only marks flagged tips: augur's comment at that call says the clock filter "will mark, but not remove bad tips", and v0 `TreeTime.clock_filter()` sets `bad_branch` and keeps the tip.

v1 follows v0: the clock filter of `timetree` flags tips (`fn filter_clock_outliers()` [packages/treetime/src/timetree/pre_loop.rs#L166](../../packages/treetime/src/timetree/pre_loop.rs#L166)), the flagged tips become bad branches (`fn bad_leaves()` [packages/treetime/src/timetree/inference/bad_branches.rs#L7](../../packages/treetime/src/timetree/inference/bad_branches.rs#L7)), and they stay in every output tree. A run that reproduces a Nextstrain build keeps tips that augur removes; ncov, dengue, rsv and seasonal-flu pass `--clock-filter-iqd` (8 for ncov, 4 for the others). On `data/sc2/4500` the clock filter flags at least 270 tips.

The maintainer approved adding the pruning as an option: off by default, so the default run keeps the v0 behavior, with a list of tips that stay even when flagged.

## Fix direction

- A setting of `timetree` that removes the flagged tips after the clock filter and before the first time inference
- A setting with the names of tips that stay in the tree even when flagged, as names and as a file with one name per line, like `--keep-ids`
- A run warning that names the removed tips
- Candidates for reuse: the node removal of `treetime prune --remove` (`fn prune_nodes()` [packages/treetime/src/prune/prune.rs#L13](../../packages/treetime/src/prune/prune.rs#L13)), the name-list reader of `--remove-file` (`name_list` in `packages/treetime-io/src/`) and `fn key_by_name()` for the names the user types

> [!IMPORTANT]
> **Decision required.** The option is approved; these parts of its design are open:
>
> - The names of the two settings, for example `--prune-outliers` and `--keep-tips` with `--keep-tips-file`, or augur's `--keep-ids`
> - Whether `clock` gets the same option, since it writes trees too, or only `timetree`
> - How removed tips appear in the outputs: absent from every tree and table, or listed in the outlier output with a "removed" flag
> - Where the removal sits in the timetree schedule: augur removes the tips, then TreeTime reroots and fits again with the tips gone; in v1 the reroot after the filter already leaves the outliers out

> [!IMPORTANT]
> **Investigation required.** `fn prune_nodes()` works on sparse reconstructions. Check whether the timetree state at the clock filter (the branch model and the partitions) can reuse it, or which part of it must move to a shared place, and how the removal of a tip collapses its parent when one child remains.

## Validation

- With the option off, the outputs of every smoke case are unchanged
- With the option on, the output tree has the tips of the input tree minus the flagged tips that the keep list does not name
- On a dataset with a clock filter, the removed tips equal the tips that `augur refine --clock-filter-iqd` with the same threshold removes, when both flag the same outliers
- A tip on the keep list stays in the tree and keeps its flag
