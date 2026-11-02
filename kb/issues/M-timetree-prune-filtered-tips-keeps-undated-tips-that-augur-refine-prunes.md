# `--prune-filtered-tips` keeps the undated tips that augur refine prunes

## augur and v0 behavior

- **v0 marks undated tips as bad**: `ClockTree` sets `bad_branch = True` on every terminal node without a date constraint and on every node whose date it cannot read ([packages/legacy/treetime/treetime/clock_tree.py#L120-L148](../../packages/legacy/treetime/treetime/clock_tree.py#L120-L148))
- **augur prunes every bad tip**: with `--clock-filter-iqd`, `def refine()` runs the clock filter and then prunes each leaf with `bad_branch` set that `--keep-ids` does not name [[src](https://github.com/nextstrain/augur/blob/58f329042e7fc631ca275cf2df44c8720f80949f/augur/refine.py#L52-L62)], so the undated tips leave the tree together with the clock outliers. Without `--clock-filter-iqd`, augur prunes nothing, and v0 dates the undated tips from the tree

## v1 behavior

`--prune-filtered-tips` removes the tips that the clock filter flags. A tip without a usable date is not flagged; the time inference dates it from the tree ([kb/decisions/timetree-uncertain-leaf-dates-are-inferred.md](../decisions/timetree-uncertain-leaf-dates-are-inferred.md)).

`--prune-undated-tips` removes the tips without a usable date: a tip without a date in the metadata, or with a value that is not a date, a numeric date, or a date range. `timetree` removes them before the time inference, together with the flagged tips when `--prune-filtered-tips` is also set; `clock` removes them before the final clock fit. The setting does not require the clock filter. The tips that `--keep-tips` and `--keep-tips-file` name stay in the tree, as with augur `--keep-ids`, and a run warning (`undated-tips-removed`) names the removed tips. The core is `PruneTips` in [packages/treetime/src/clock/prune_tips.rs](../../packages/treetime/src/clock/prune_tips.rs).

An augur workflow with `--clock-filter-iqd` matches augur refine when it passes both settings. A six-tip tree whose tip `F` has no date keeps `F` under `--clock-filter 4 --prune-filtered-tips`, and loses it when `--prune-undated-tips` is added:

```bash
printf '((A:0.002,B:0.003):0.001,((C:0.002,D:0.004):0.002,(E:0.003,F:0.001):0.001):0.001);\n' > tree.nwk
printf 'strain\tdate\nA\t2020-01-01\nB\t2021-01-01\nC\t2020-06-01\nD\t2022-01-01\nE\t2021-06-01\nF\t\n' > metadata.tsv
treetime timetree --tree tree.nwk --metadata metadata.tsv --sequence-length 1000 --clock-rate 0.001 \
  --clock-filter 4 --prune-filtered-tips --prune-undated-tips --keep-root --branch-length-mode input --seed 1 \
  --output-tree-nwk out.nwk --output-node-data out.json
```

On replayed Nextstrain builds that use `--clock-filter-iqd`, the trees of v0 and v1 keep different tips, for example 12 tips that only v1 keeps in the seasonal-cov 229E build and 5 in the rubella genome build. Steps after `refine` that read dates from the metadata, such as `augur frequencies`, can then meet tips without a date.

> [!IMPORTANT]
> **Investigation required.** How many of the tips that only v1 keeps in these builds lack a date, and how many come from a different outlier set of the clock filter ([M-clock-filter-residual-parity.md](M-clock-filter-residual-parity.md)); and whether `augur frequencies` stops on the kept undated tips of the rubella and seasonal-cov workflows.
