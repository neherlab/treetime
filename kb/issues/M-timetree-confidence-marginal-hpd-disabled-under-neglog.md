# Timetree confidence intervals drop the marginal-posterior HPD contribution under NegLog

Node confidence intervals combine two sources: the highest-posterior-density (HPD) region of the marginal time distribution (mutation stochasticity) and the rate-susceptibility interval (clock-rate uncertainty). The marginal-posterior source is disabled: `extract_confidence_intervals` in `packages/treetime/src/timetree/confidence.rs` hardcodes `mutation_contribution = None` and comments out the `hpd_region` call, because `hpd_region` reads plain-probability ordinates and the time distribution now stores negative-log ordinates (`Distribution<NegLog>`).

## Impact and scope

- A node whose only confidence source is its marginal time distribution (no rate-susceptibility dates) falls back to the point estimate, so its interval collapses to `[date, date]`.
- A node with both sources reports only the rate contribution, so the combined interval is never wider than the rate source alone.
- v0 uses `get_max_posterior_region(fraction=0.9)`, the narrowest interval holding 90% of the mass, which is narrower than an equal-tailed interval for skewed posteriors. v1 currently emits no marginal contribution at all.

## Root cause

`hpd_region` integrates a probability density and is defined for `Distribution<Plain>` only. Under `NegLog` the stored ordinate is `-ln p`, so the plain-space integration does not apply. A NegLog-aware HPD (integrate after converting through the peak-normalized `to_plain_normalized`, or integrate directly in neg-log space) has not been implemented.

## Tests

Three tests in `packages/treetime/src/commands/timetree/output/__tests__/test_confidence_extract.rs` exercise the marginal-posterior HPD and are marked `#[ignore = "marginal-posterior HPD disabled pending NegLog-aware HPD"]`:

- `test_extract_confidence_intervals_with_distribution`
- `test_extract_confidence_intervals_combined_wider_than_either`
- `test_extract_confidence_intervals_skewed_distribution_hpd`

The ignored tests already store their distributions on the neg-log axis, so they are ready to re-enable once the HPD path returns.

## Fix approach

Implement a NegLog-aware HPD region and restore the `mutation_contribution` branch in `extract_confidence_intervals`, then remove the `#[ignore]` from the three tests. Tracked under the log-space distribution work in [kb/proposals/distribution-log-space-and-hard-soft-boundaries.md](../proposals/distribution-log-space-and-hard-soft-boundaries.md) (HPD is Part B/D there).
