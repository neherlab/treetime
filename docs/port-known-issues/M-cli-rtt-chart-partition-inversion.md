# RTT chart gather_points swaps normal and outlier series

`gather_points` at [packages/treetime/src/cli/rtt_chart.rs#L133-L177](../../packages/treetime/src/cli/rtt_chart.rs#L133-L177) calls `results.iter().partition(|r| r.is_outlier)` on line 139. Rust's `Iterator::partition` returns `(matches_predicate, rest)`, so the binding `let (norms, outliers) = ...` assigns outlier results to `norms` and non-outlier results to `outliers`.

Lines 141-150 then build `norm_points` from `outliers.into_iter()` and `outlier_points` from `norms.into_iter()`, creating a second inversion that cancels the first for the populated `PointsResult` fields. The double-swap means the final chart output is accidentally correct today.

## Latent defects

1. The intermediate variable names (`norms`, `outliers`) are semantically inverted, so any refactoring that uses them at face value produces wrong-color legends and reversed series.

2. The first iterator (line 141) uses `.cloned()` while the second (line 147) does not, creating an asymmetry in copy semantics between the two branches. The first branch clones `ClockRegressionResult` structs; the second moves references. Both happen to work because `filter_map` extracts `(f32, f32)` tuples downstream, but the inconsistency is a maintenance trap.

3. The function uses `assert!(!results.is_empty())` (line 137) and `.minmax().into_option().unwrap()` (lines 154, 165), violating the project's "NEVER panic" rule. Empty `results` or all-`None` dates cause panics.

## Affected code

- Partition binding: [packages/treetime/src/cli/rtt_chart.rs#L139](../../packages/treetime/src/cli/rtt_chart.rs#L139)
- First iterator (with `.cloned()`): [packages/treetime/src/cli/rtt_chart.rs#L141-L145](../../packages/treetime/src/cli/rtt_chart.rs#L141-L145)
- Second iterator (without `.cloned()`): [packages/treetime/src/cli/rtt_chart.rs#L147-L150](../../packages/treetime/src/cli/rtt_chart.rs#L147-L150)
- Panic sites: [packages/treetime/src/cli/rtt_chart.rs#L137](../../packages/treetime/src/cli/rtt_chart.rs#L137), [rtt_chart.rs#L154](../../packages/treetime/src/cli/rtt_chart.rs#L154), [rtt_chart.rs#L165](../../packages/treetime/src/cli/rtt_chart.rs#L165)

## Fix

Swap the partition binding to `let (outliers, norms) = ...`, remove the compensating swap in the lambda builders, unify `.cloned()` usage on both branches, replace `assert!` with `make_error!`, and replace `.unwrap()` with `.ok_or_else(|| make_report!(...))`.

Add a unit test for `gather_points` that asserts a single outlier appears in `outlier_points` (not in `norm_points`), a single normal appears in `norm_points`, and that empty input returns an error.
