# Timetree confidence intervals drop the marginal-posterior HPD contribution under NegLog

Node confidence intervals combine two sources: the highest-posterior-density (HPD) region of the marginal time distribution (mutation stochasticity) and the rate-susceptibility interval (clock-rate uncertainty). The marginal-posterior source is missing: `extract_confidence_intervals` sets `mutation_contribution` to `None` [`packages/treetime/src/timetree/confidence.rs#L137`](../../packages/treetime/src/timetree/confidence.rs#L137), because v1 has no HPD region for the negative-log ordinates the time distribution stores (`Distribution<NegLog>`).

## Impact and scope

- A node whose only confidence source is its marginal time distribution (no rate-susceptibility dates) falls back to the point estimate, so its interval collapses to `[date, date]`.
- A node with both sources reports only the rate contribution, so the combined interval is never wider than the rate source alone.
- v0 uses `get_max_posterior_region(fraction=0.9)`, the narrowest interval holding 90% of the mass, which is narrower than an equal-tailed interval for skewed posteriors. v1 currently emits no marginal contribution at all.

## Root cause

An HPD region integrates a probability density. Under `NegLog` the stored ordinate is `-ln p`, so a region must either convert to peak-normalized plain probabilities before integrating or integrate directly in neg-log space. Neither is implemented.

## Tests

Three tests in `packages/treetime/src/timetree/__tests__/test_confidence_extract.rs` exercise the marginal-posterior HPD and are marked `#[ignore = "marginal-posterior HPD disabled pending NegLog-aware HPD"]`:

- `test_extract_confidence_intervals_with_distribution`
- `test_extract_confidence_intervals_combined_wider_than_either`
- `test_extract_confidence_intervals_skewed_distribution_hpd`

The ignored tests already store their distributions on the neg-log axis, so they are ready to re-enable once the HPD path returns.

## Fix approach

Implement a NegLog-aware HPD region and pass its interval as `mutation_contribution` in `extract_confidence_intervals`, then remove the `#[ignore]` from the three tests. Tracked under the log-space distribution work in [kb/proposals/distribution-log-space-and-hard-soft-boundaries.md](../proposals/distribution-log-space-and-hard-soft-boundaries.md) (HPD is Part B/D there).
