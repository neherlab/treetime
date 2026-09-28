# Fix crash and extra leaf when rerooting a tree whose root has one child

Source: [neherlab/treetime#959](https://github.com/neherlab/treetime/issues/959). Analysis, reproduction data, and rejected options: [kb/reports/reroot-single-child-root.md](../reports/reroot-single-child-root.md).

## Problem

Rerooting a tree whose root has one child turns the old root into an undated leaf. `treetime` and `treetime clock` then crash in the root-to-tip plot, because the regression's tip value and the plot's outlier mask test the `numpy` flag `bad_branch` in different ways. Without the crash, the extra leaf would be written to the output trees.

`resolve_polytomies()` also removes one-child nodes that carry a date constraint, so their dates are dropped without a message.

## Work items

- Remove an undated one-child root in `TreeTime.reroot()` before rerooting. Keep a dated one-child root. Do not change `--keep-root` runs
- Store `bad_branch` as a Python `bool` in `ClockTree._assign_dates()` and test it by truthiness in the regression tip value, `init_date_constraints()`, and `convert_dates()`
- After a reroot, raise an error if a leaf without a date was created
- Keep one-child nodes with a date constraint in `resolve_polytomies()`
- Add tests without network data for each item, and for the `timetree` and `clock` command paths of the reported case

## Done when

- The commands in the report exit with 0, and the output trees contain exactly the input samples
- `dates.tsv` for the reported tree equals the one for the same tree without the one-child root
- `--keep-root` output is unchanged
- A dated one-child node inside the tree keeps its date in the default run
