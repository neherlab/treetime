# `--prune-filtered-tips` keeps the undated tips that augur refine prunes

> [!IMPORTANT]
> **Decision required.** augur refine with a clock filter removes every tip that v0 marks as bad, which includes the tips without a usable date. v1 `--prune-filtered-tips` removes only the tips that the clock filter flags, and keeps the undated tips with an inferred date. See "Decision required" below.

## augur and v0 behavior

- **v0 marks undated tips as bad**: `ClockTree` sets `bad_branch = True` on every terminal node without a date constraint and on every node whose date it cannot read ([packages/legacy/treetime/treetime/clock_tree.py#L120-L148](../../packages/legacy/treetime/treetime/clock_tree.py#L120-L148))
- **augur prunes every bad tip**: with `--clock-filter-iqd`, `def refine()` runs the clock filter and then prunes each leaf with `bad_branch` set that `--keep-ids` does not name [[src](https://github.com/nextstrain/augur/blob/58f329042e7fc631ca275cf2df44c8720f80949f/augur/refine.py#L52-L62)], so the undated tips leave the tree together with the clock outliers. Without `--clock-filter-iqd`, augur prunes nothing, and v0 dates the undated tips from the tree

## v1 behavior

`fn prune_filtered_tips()` [packages/treetime/src/timetree/pipeline.rs#L189-L229](../../packages/treetime/src/timetree/pipeline.rs#L189-L229) removes the tips that the clock filter flags. A tip without a date is not flagged; the time inference dates it from the tree ([kb/decisions/timetree-uncertain-leaf-dates-are-inferred.md](../decisions/timetree-uncertain-leaf-dates-are-inferred.md)). The help text of the flag says that the run removes the flagged tips "as augur refine does" ([packages/app-commands/src/commands/shared/clock_filter.rs#L52](../../packages/app-commands/src/commands/shared/clock_filter.rs#L52)), which is true for the clock outliers only.

A six-tip tree whose tip `F` has no date keeps `F`, with the inferred date 2019-08-31, under `--clock-filter 4 --prune-filtered-tips`:

```bash
printf '((A:0.002,B:0.003):0.001,((C:0.002,D:0.004):0.002,(E:0.003,F:0.001):0.001):0.001);\n' > tree.nwk
printf 'strain\tdate\nA\t2020-01-01\nB\t2021-01-01\nC\t2020-06-01\nD\t2022-01-01\nE\t2021-06-01\nF\t\n' > metadata.tsv
treetime timetree --tree tree.nwk --metadata metadata.tsv --sequence-length 1000 --clock-rate 0.001 \
  --clock-filter 4 --prune-filtered-tips --keep-root --branch-length-mode input --seed 1 \
  --output-tree-nwk out.nwk --output-node-data out.json
```

On replayed Nextstrain builds that use `--clock-filter-iqd`, the trees of v0 and v1 keep different tips, for example 12 tips that only v1 keeps in the seasonal-cov 229E build and 5 in the rubella genome build. Steps after `refine` that read dates from the metadata, such as `augur frequencies`, can then meet tips without a date.

> [!IMPORTANT]
> **Investigation required.** How many of the tips that only v1 keeps in these builds lack a date, and how many come from a different outlier set of the clock filter ([M-clock-filter-residual-parity.md](M-clock-filter-residual-parity.md)); and whether `augur frequencies` stops on the kept undated tips of the rubella and seasonal-cov workflows.

## Decision required

- **Prune the undated tips with the flagged ones**: `--prune-filtered-tips` removes every tip without a usable date as well, as augur refine with a clock filter does. Changes which tips the outputs contain
- **A separate setting**: keep `--prune-filtered-tips` as it is, and add a setting that removes undated tips, which augur workflows then pass. Two settings for one augur behavior
- **Keep the behavior and correct the help text**: undated tips stay in every run; augur workflows that need them out must filter the metadata before `refine`
