# Clock regression keeps the dates of clock-filter outliers

The clock filter marks leaves whose root-to-tip residual exceeds the IQD threshold as outliers. The regression is meant to leave these leaves out, as v0 does: v0 gives bad-branch tips no value in the tree regression (`tip_value` returns `None` for them, [packages/legacy/treetime/treetime/clock_tree.py#L275](../../packages/legacy/treetime/treetime/clock_tree.py#L275)).

In v1 the outlier set only changes the message a leaf stores for its parent edge (`clock_to_parent`). The message the parent sums (`clock_from_child`) is built from the leaf date regardless of the outlier set ([packages/treetime/src/clock/clock_regression.rs#L281](../../packages/treetime/src/clock/clock_regression.rs#L281), [packages/treetime/src/clock/clock_regression.rs#L305-L307](../../packages/treetime/src/clock/clock_regression.rs#L305-L307)). The root statistics therefore include every outlier date, and the outlier flag only affects root candidates on the outlier's own branch during reroot.

## Evidence

`test_clock_regression_outliers_are_left_out_of_the_fit` in `packages/treetime/src/clock/__tests__/test_clock_regression.rs` (ignored) fits a five-leaf tree twice with the root kept: once with leaf E marked as outlier, once with E undated. The two clock rates differ (`3.506e-3` with E as outlier, `7.373e-3` with E undated), while v0 gives both the same rate because it drops the date in both cases.

## Impact

`--clock-filter` in `treetime clock` and `treetime timetree` flags outliers and reports them, but the clock rate, intercept and R² still use their dates. This affects every clock model fitted after the filter, the timetree rounds included. It may explain part of the R² difference on dengue/100 noted in [M-clock-filter-residual-parity.md](M-clock-filter-residual-parity.md).

## Open question

Fixing the backward message changes the clock output of every run with flagged outliers, so it is an output change that needs a parity check against v0 on the smoke datasets.
