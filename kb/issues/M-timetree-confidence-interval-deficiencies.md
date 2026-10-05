# Timetree confidence interval computation deficiencies

> [!IMPORTANT]
> **Decision required.** Of the three defects this issue was opened for, one no longer exists, one matches v0, and one is a v1-only step that differs from its first description:
>
> - Rate susceptibility no longer mutates the graph: `fn compute_rate_susceptibility()` [`packages/treetime/src/timetree/confidence.rs#L23-L72`](../../packages/treetime/src/timetree/confidence.rs#L23-L72) runs `run_timetree` three times on copies of `TimeInferenceInputs` with scaled gammas; the inputs hold the graph as `&Graph` ([`packages/treetime/src/timetree/inference/runner.rs#L79-L80`](../../packages/treetime/src/timetree/inference/runner.rs#L79-L80)). No restoration pass exists to verify
> - The symmetric rate term matches v0: v0 computes `c + x * np.abs(y - c)` ([`packages/legacy/treetime/treetime/clock_tree.py#L1085`](../../packages/legacy/treetime/treetime/clock_tree.py#L1085)), the same as v1 [`packages/treetime/src/timetree/confidence.rs#L129-L130`](../../packages/treetime/src/timetree/confidence.rs#L129-L130). Under the reference-parity rule this is expected behavior unless a v0 erratum is recorded
> - The last step widens the interval to contain the point estimate (`lower.min(date)`, `upper.max(date)`, [`packages/treetime/src/timetree/confidence.rs#L110-L111`](../../packages/treetime/src/timetree/confidence.rs#L110-L111)); it does not move the point estimate. v0 has no such step: `get_confidence_interval` returns the `combine_confidence` result directly ([`packages/legacy/treetime/treetime/clock_tree.py#L1128-L1145`](../../packages/legacy/treetime/treetime/clock_tree.py#L1128-L1145))
>
> Options: close the issue; keep it only for the v1-only widening step (remove the step for parity, or record it as an intentional change); or keep the asymmetry question open as a proposed v0 erratum, which needs evidence.

## Summary

Confidence intervals of `timetree --confidence` combine a rate-susceptibility term with the bounds of the node time distribution. The computation differs from v0 and from augur in the places listed below.

## Details

### Symmetric rate contribution

`fn date_uncertainty_due_to_rate()` [`packages/treetime/src/timetree/confidence.rs#L125-L132`](../../packages/treetime/src/timetree/confidence.rs#L125-L132)

Both bounds are computed from the absolute difference between the rate-perturbed date and the nominal date. This forces the rate contribution to be symmetric around the point estimate, regardless of whether the rate-date relationship is asymmetric (as it can be for short branches where the Poisson likelihood is skewed). v0 does the same ([`packages/legacy/treetime/treetime/clock_tree.py#L1085`](../../packages/legacy/treetime/treetime/clock_tree.py#L1085)).

### Interval widened to contain the point estimate

`fn extract_confidence_intervals()` [`packages/treetime/src/timetree/confidence.rs#L110-L111`](../../packages/treetime/src/timetree/confidence.rs#L110-L111)

After combining the contributions, the code widens `[lower, upper]` so that it contains the reported date. A case where the computed interval and the reported date are inconsistent (e.g. from the symmetric approximation above, or numerical issues in rate susceptibility) is hidden instead of reported. v0 has no equivalent step.

### Mutation contribution is not used

`fn extract_confidence_intervals()` fixes `mutation_contribution` to `None` ([`packages/treetime/src/timetree/confidence.rs#L93`](../../packages/treetime/src/timetree/confidence.rs#L93)). v0 passes the marginal inverse-CDF quantiles of the node time as `c2` to `combine_confidence` ([`packages/legacy/treetime/treetime/clock_tree.py#L1138`](../../packages/legacy/treetime/treetime/clock_tree.py#L1138)). Without a rate contribution, v1 reports the degenerate interval `[date, date]`.

### num_date_confidence composition diverges from augur's marginal HPD

`num_date_confidence` in `timetree.augur-node-data.json` (and `auspice_tree.json` `num_date.confidence`) is the `[lower, upper]` 90% region produced by this module. Augur sets `num_date_confidence = list(tt.get_max_posterior_region(n, 0.9))`, the pure marginal-posterior HPD. v1 instead reports the rate-susceptibility contribution, clamped to the support of the node time distribution and made symmetric by the `abs()` above. The field mapping is correct (a 90% region), but the bounds can differ numerically from augur. This is a consumed field (auspice colors and HPD bars use it), so closing this affects auspice output, not just the node data file.

## Impact

- Rate contributions are symmetric even where the rate-date relationship is asymmetric
- Inconsistencies between the interval and the point estimate are hidden instead of flagged
- `num_date_confidence` differs from augur's marginal HPD
