# Clock pre-filter allows negative rate during root finding

## Deviation

v1's clock pre-filter step uses `force_positive_rate: false` when finding the best root for outlier detection. v0 uses `force_positive=True` (the default) for the same step.

## Rationale

v1's clock regression produces negative estimated rates at all root positions for some datasets (e.g. dengue/100 with 8-10 outliers) where v0 finds at least one positive-rate root. The root cause of this difference in regression results is not fully understood (tracked separately).

Without this relaxation, v1 crashes with "Clock rate is negative for all root positions" before outlier filtering runs. The pre-filter clock model only needs to be good enough for IQD-based outlier detection, which uses absolute deviations and is slope-sign-invariant for extreme outliers.

After filtering, the final clock model estimation warns on a non-positive rate but continues, matching v0. The `--allow-negative-rate` flag relaxes the root search to accept negative-slope roots. See [timetree-rejects-negative-clock-rate.md](timetree-rejects-negative-clock-rate.md) for the contrasting timetree behavior.

## Impact

- Pre-filter may select a different root position than v0 for the initial outlier detection pass
- Outlier sets may differ between v0 and v1 (additional differences from `clock_filter_inplace` vs v0's `residual_filter` are tracked separately in known issues)
- Final clock model quality depends on the effectiveness of the pre-filter, but the alternative (crashing) is strictly worse

## v0 Reference

`packages/legacy/treetime/treetime/treetime.py:486-489`: `self.reroot(root='least-squares', covariation=False)` uses default `force_positive=True`.

`packages/legacy/treetime/treetime/treeregression.py:321-353`: `find_best_root(force_positive=True)` only accepts positive-slope roots.

`packages/legacy/treetime/treetime/wrappers.py:984`: `myTree.reroot(params.reroot, force_positive=not params.allow_negative_rate)` applies the user flag to the FINAL reroot, not the pre-filter.

## v1 Implementation

`packages/treetime/src/commands/clock/run.rs`: `estimate_clock_model_with_prefilter()` uses `RerootParams { force_positive_rate: false }` for the pre-filter step. The final step uses `force_positive_rate: !allow_negative_rate`.
