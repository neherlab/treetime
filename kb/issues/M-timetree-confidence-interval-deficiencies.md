# Timetree confidence interval computation deficiencies

> [!IMPORTANT]
> **Decision required.** Of the three defects this issue was opened for, one no longer exists, one matches v0, and one is a v1-only step that differs from its first description:
>
> - Rate susceptibility no longer mutates the graph: `fn compute_rate_susceptibility()` [`packages/treetime/src/timetree/confidence.rs#L25-L74`](../../packages/treetime/src/timetree/confidence.rs#L25-L74) runs `run_timetree` three times on copies of `TimeInferenceInputs` with scaled gammas; the inputs hold the graph as `&Graph` ([`packages/treetime/src/timetree/inference/runner.rs#L79-L80`](../../packages/treetime/src/timetree/inference/runner.rs#L79-L80)). No restoration pass exists to verify
> - The symmetric rate term matches v0: v0 computes `c + x * np.abs(y - c)` ([`packages/legacy/treetime/treetime/clock_tree.py#L1085`](../../packages/legacy/treetime/treetime/clock_tree.py#L1085)), the same as v1 [`packages/treetime/src/timetree/confidence.rs#L166-L167`](../../packages/treetime/src/timetree/confidence.rs#L166-L167). Under the reference-parity rule this is expected behavior unless a v0 erratum is recorded
> - The last step widens the interval to contain the point estimate (`lower.min(date)`, `upper.max(date)`, [`packages/treetime/src/timetree/confidence.rs#L119-L120`](../../packages/treetime/src/timetree/confidence.rs#L119-L120)); it does not move the point estimate. v0 has no such step: `get_confidence_interval` returns the `combine_confidence` result directly ([`packages/legacy/treetime/treetime/clock_tree.py#L1128-L1145`](../../packages/legacy/treetime/treetime/clock_tree.py#L1128-L1145))
>
> Options: close the issue; keep it only for the v1-only widening step (remove the step for parity, or record it as an intentional change); or keep the asymmetry question open as a proposed v0 erratum, which needs evidence.

## Summary

Confidence intervals of `timetree --confidence` combine a rate-susceptibility term with the bounds of the node time distribution. The computation differs from v0 and from augur in the places listed below.

## Details

### Symmetric rate contribution

`fn date_uncertainty_due_to_rate()` [`packages/treetime/src/timetree/confidence.rs#L162-L169`](../../packages/treetime/src/timetree/confidence.rs#L162-L169)

Both bounds are computed from the absolute difference between the rate-perturbed date and the nominal date. This forces the rate contribution to be symmetric around the point estimate, regardless of whether the rate-date relationship is asymmetric (as it can be for short branches where the Poisson likelihood is skewed). v0 does the same ([`packages/legacy/treetime/treetime/clock_tree.py#L1085`](../../packages/legacy/treetime/treetime/clock_tree.py#L1085)).

### Interval widened to contain the point estimate

`fn node_confidence_interval()` [`packages/treetime/src/timetree/confidence.rs#L119-L120`](../../packages/treetime/src/timetree/confidence.rs#L119-L120)

After combining the contributions, the code widens `[lower, upper]` so that it contains the reported date. A case where the computed interval and the reported date are inconsistent (e.g. from the symmetric approximation above, or numerical issues in rate susceptibility) is hidden instead of reported. v0 has no equivalent step.

### num_date_confidence differs from augur only by the widening step

`num_date_confidence` in `timetree.augur-node-data.json` (and `auspice_tree.json` `num_date.confidence`) is the `[lower, upper]` 90% region produced by this module. Augur sets `num_date_confidence = list(tt.get_max_posterior_region(n, 0.9))`. v1 composes the region the same way: the maximum posterior region of the marginal time posterior, combined in quadrature with the rate-susceptibility term when one exists, clipped to the support of the posterior. A node with a point posterior (for example a tip with an exact date) gets `[date, date]` with its reported date, as augur does through v0 `get_max_posterior_region` ([`packages/legacy/treetime/treetime/clock_tree.py#L1170-L1171`](../../packages/legacy/treetime/treetime/clock_tree.py#L1170-L1171)). The bounds differ from augur only where the widening step above moves a bound to the point estimate, and on the posteriors of [kb/v0-errata/timetree-max-posterior-region-multimodal-side.md](../v0-errata/timetree-max-posterior-region-multimodal-side.md) and [kb/v0-errata/timetree-max-posterior-region-brent-plateau.md](../v0-errata/timetree-max-posterior-region-brent-plateau.md). This is a consumed field (auspice colors and HPD bars use it), so a change of the widening step affects auspice output, not just the node data file.

## Impact

- Rate contributions are symmetric even where the rate-date relationship is asymmetric
- Inconsistencies between the interval and the point estimate are hidden instead of flagged
- `num_date_confidence` differs from augur where the widening step applies
