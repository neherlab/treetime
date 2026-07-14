# No v0 golden master test for min-dev reroot

> [!IMPORTANT]
> **Investigation required.** End-to-end parity remains unverified until an independent v0 golden master covers edge selection and split position.

There is no test comparing the v1 min-dev root position against a v0 oracle. The min-dev objective (minimize root-to-tip divergence variance) exists in v0, so a golden master is feasible, but none is captured.

## Background

v0 implements min-dev rerooting in `def TreeRegression.optimal_reroot()` [packages/legacy/treetime/treetime/treeregression.py](../../packages/legacy/treetime/treetime/treeregression.py) with `slope=0`, which minimizes the variance of root-to-tip distances. v1 reproduces this with `struct DivStats` divergence-only scoring in the generic reroot module. [kb/proposals/optimize-reroot-support.md](../proposals/optimize-reroot-support.md) lists as a validation criterion that the root from `optimize --reroot=min-dev` matches the root from v0 min-dev rerooting on the same tree.

Current coverage validates the scoring algebra in isolation: `struct DivStats` is cross-checked against `fn ClockSet::propagate_averages()` in unit tests [packages/treetime/src/reroot/__tests__/test_div_stats.rs](../../packages/treetime/src/reroot/__tests__/test_div_stats.rs). No end-to-end test confirms that the selected root edge and split position match v0 on a real dataset.

## Impact

Negligible. The objective is unit-tested against the clock-set propagation it mirrors, so the risk is a wiring-level discrepancy (edge selection, endpoint handling, branch-length conventions) rather than a scoring error. A golden master would confirm end-to-end v0 parity.

## Proposed test

Capture the v0 min-dev root on a small dataset by running v0 rerooting, then assert that the v1 optimize path selects the same root edge and matches the split fraction and incident lengths within the project's $10^{-6}$ numerical contract. Compare the clock path separately as an expected detector for its known objective and split-optimizer divergences.

## Locations

- Generic reroot implementation [packages/treetime/src/reroot/mod.rs](../../packages/treetime/src/reroot/mod.rs)
- Existing algebra cross-check [packages/treetime/src/reroot/__tests__/test_div_stats.rs](../../packages/treetime/src/reroot/__tests__/test_div_stats.rs)
- v0 reference [packages/legacy/treetime/treetime/treeregression.py](../../packages/legacy/treetime/treetime/treeregression.py)
- [kb/proposals/reroot-generic-scoring-architecture.md](../proposals/reroot-generic-scoring-architecture.md)
- [kb/proposals/optimize-reroot-support.md](../proposals/optimize-reroot-support.md)

## Related tickets

- [kb/tickets/test-reroot-min-dev-v0-golden-master.md](../tickets/test-reroot-min-dev-v0-golden-master.md)

## Related issues

- [M-clock-mindev-wrong-objective.md](M-clock-mindev-wrong-objective.md)
- [N-reroot-split-optimizer-default-diverges-from-v0.md](N-reroot-split-optimizer-default-diverges-from-v0.md)
