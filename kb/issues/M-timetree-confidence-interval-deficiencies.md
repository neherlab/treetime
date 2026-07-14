# Timetree confidence interval computation deficiencies

## Summary

Three defects in the confidence interval computation: the rate susceptibility function mutates graph state without postcondition assertion, date uncertainty is forced symmetric, and the CI clamp hides inconsistencies.

## Details

### compute_rate_susceptibility mutates graph 3 times without postcondition

`packages/treetime/src/timetree/confidence.rs:56-124:`

The function runs marginal inference three times (upper rate, lower rate, restored original rate), relying on the third run to restore the graph to its pre-call state. No postcondition assertion verifies that restoration succeeded. If the third pass fails or produces different rounding, subsequent operations read corrupted state.

### date_uncertainty_due_to_rate uses abs() making CI symmetric

`packages/treetime/src/timetree/confidence.rs:155-156:`

Both upper and lower confidence bounds are computed as the absolute difference between the rate-perturbed date and the nominal date. This forces the confidence interval to be symmetric around the point estimate, regardless of whether the rate-date relationship is asymmetric (which it is for short branches where the Poisson likelihood is skewed).

### CI clamp forces point estimate inside interval

`packages/treetime/src/timetree/confidence.rs:219-220:`

After computing confidence bounds, the code clamps the point estimate to lie within `[lower, upper]`. This hides cases where the CI computation and the reported date are inconsistent (e.g., due to the symmetric approximation above or numerical issues in rate susceptibility). The inconsistency should be reported, not masked.

### num_date_confidence composition diverges from augur's marginal HPD

`num_date_confidence` in `timetree.augur-node-data.json` (and `auspice_tree.json` `num_date.confidence`) is the `[lower, upper]` 90% region produced by this module. Augur sets `num_date_confidence = list(tt.get_max_posterior_region(n, 0.9))`, the pure marginal-posterior HPD. v1 combines the marginal HPD with a rate-susceptibility contribution (quadrature) when `--confidence` runs with rate uncertainty, and forces symmetry via the `abs()` above. The field mapping is correct (a 90% region), but the bounds can differ numerically from augur. This is a consumed field (auspice colors and HPD bars use it), so closing this affects auspice output, not just the node data file.

## Impact

- Confidence intervals are symmetric when they should be asymmetric for short branches
- Graph state corruption if third rate-susceptibility pass diverges from original
- Inconsistent point estimates hidden instead of flagged

## Related tickets

- [kb/tickets/timetree-confidence-interval-computation-deficiencies.md](../tickets/timetree-confidence-interval-computation-deficiencies.md)
