# Add an independent end-to-end oracle for minimum-deviation rerooting

Validate generic minimum-deviation rerooting against deterministic analytical fixtures rather than TreeTime v0's defective candidate score.

## Acceptance criteria

- Fixtures have independently derived minimum-variance root edges, split fractions, and resulting incident branch lengths.
- The generic reroot path matches each analytical oracle with absolute error at most $10^{-6}$.
- The clock path selects the same edge and split under the same no-covariation weighting, where normalized and unnormalized fixed-zero-rate scores differ only by a root-independent factor.
- Endpoint roots and ties have explicit deterministic expectations.
- Expected values are derived analytically and are not computed by v1 helpers.
- Perturbing an expected edge or split fraction makes the test fail clearly.
- Any v0 comparison is labeled as regression coverage for [kb/v0-errata/clock-min-dev-fixed-slope-score.md](../v0-errata/clock-min-dev-fixed-slope-score.md), not as a correctness oracle.

## Related issues

- Source: [kb/issues/N-reroot-missing-min-dev-end-to-end-oracle.md](../issues/N-reroot-missing-min-dev-end-to-end-oracle.md)
- [kb/issues/N-reroot-split-optimizer-default-diverges-from-v0.md](../issues/N-reroot-split-optimizer-default-diverges-from-v0.md)
