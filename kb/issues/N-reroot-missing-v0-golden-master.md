# No v0 golden master test for min-dev reroot

There is no test comparing the v1 min-dev root position against a v0 oracle. The min-dev objective (minimize root-to-tip divergence variance) exists in v0, so a golden master is feasible, but none is captured.

## Background

v0 implements min-dev rerooting in `TreeRegression.optimal_reroot` with `slope=0`, which minimizes the variance of root-to-tip distances (`packages/legacy/treetime/treetime/treeregression.py`, `_optimal_root_along_branch`). v1 reproduces this with `DivStats` divergence-only scoring in the generic reroot module. The optimize-reroot proposal lists as a validation criterion that the root from `optimize --reroot=min-dev` matches the root from v0 min-dev rerooting on the same tree.

Current coverage validates the scoring algebra in isolation: `DivStats` is cross-checked against `ClockSet::propagate_averages` in unit tests (`packages/treetime/src/reroot/__tests__/test_div_stats.rs`). No end-to-end test confirms that the selected root edge and split position match v0 on a real dataset.

## Impact

Negligible. The objective is unit-tested against the clock-set propagation it mirrors, so the risk is a wiring-level discrepancy (edge selection, endpoint handling, branch-length conventions) rather than a scoring error. A golden master would confirm end-to-end v0 parity.

## Proposed test

Capture the v0 min-dev root on a small dataset (e.g. flu/h3n2/20) by running v0 rerooting, then assert the v1 root edge and split fraction agree within tolerance. Cross-check against `clock --reroot=min-dev`, which shares the same objective.

## Locations

- `packages/treetime/src/reroot/` (min-dev scoring and search under test)
- `packages/treetime/src/reroot/__tests__/test_div_stats.rs` (existing algebra cross-check)
- v0 reference: `packages/legacy/treetime/treetime/treeregression.py`
