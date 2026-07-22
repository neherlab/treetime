# Initial coalescent Tc is a hardcoded constant instead of data-derived

> **Resolved.** `INITIAL_COALESCENT_TC` has been removed. The constant-Tc optimizer
> is now a closed-form analytic solve ($T_c = I/M$, see
> [decisions/coalescent-analytic-tc-optimization.md](../decisions/coalescent-analytic-tc-optimization.md)),
> so it needs no starting guess. If it ever fails (a degenerate tree with no
> mergers), `estimate_coalescent_tc` falls back to the previous round's $T_c$, then
> to the user-supplied `--coalescent` value, and otherwise applies no coalescent
> prior that round rather than an invented timescale. No hardcoded default remains.

`INITIAL_COALESCENT_TC` is a fixed constant (`5.0` after [PR#851](https://github.com/neherlab/treetime/pull/851), previously `0.001`) used as the fallback coalescent timescale when Brent optimization fails. The value only affects the error-recovery path -- the optimizer searches its full bracket `[-20, 2]` in log space regardless of starting point.

A data-derived initial value (e.g., root-to-tip time span from the clock regression) would adapt to dataset timescale automatically. Datasets with timescales far from 5 years (sub-annual outbreaks, deep evolutionary trees) get a poor fallback when optimization fails.

## Location

[`packages/treetime/src/timetree/pipeline.rs`](../../packages/treetime/src/timetree/pipeline.rs): `INITIAL_COALESCENT_TC` constant and its use in `CoalescentInitialization::Optimize`.

## Impact

Low. The Brent optimizer succeeds for most datasets. The fallback path matters only when the tree's coalescent structure is too degenerate for optimization to converge.
