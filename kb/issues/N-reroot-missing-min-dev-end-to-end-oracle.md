# No independent end-to-end oracle for min-dev rerooting

> [!IMPORTANT]
> **Investigation required.** Unit and clock-dispatch tests cover the scoring algebra, but the generic root search still lacks an independent end-to-end oracle for edge selection and split position.

## Background

The intended `min-dev` objective minimizes root-to-tip distance variance. `struct DivStats` implements its date-free sufficient statistics, and clock `RerootMethod::MinDev` dispatches to fixed-zero-rate WLS. Existing tests cross-check `DivStats` propagation against `ClockSet` and verify the clock dispatch against a direct `RootObjective::FixedRate(0.0)` search.

TreeTime v0 cannot serve as the correctness oracle for this objective. Although v0 passes `slope=0` and documents minimum-variance behavior, its candidate score retains the estimated-rate regression term. See [kb/v0-errata/clock-min-dev-fixed-slope-score.md](../v0-errata/clock-min-dev-fixed-slope-score.md).

## Impact

Negligible. The algebra and clock wiring are tested, leaving a residual risk in the generic search's edge selection, endpoint handling, and split-position conventions.

## Required oracle

Add deterministic analytical fixtures whose minimum-variance root edge and split can be derived independently of either v1 search implementation. Cross-check the generic and clock paths under the same no-covariation weighting, where their scores differ only by a root-independent constant factor. Keep any v0 comparison as a regression test for the documented erratum rather than as a correctness oracle.

## Acceptance criteria

- Fixtures have independently derived minimum-variance root edges, split fractions, and resulting incident branch lengths
- The generic reroot path matches each analytical oracle with absolute error at most $10^{-6}$
- The clock path selects the same edge and split under the same no-covariation weighting
- Endpoint roots and ties have explicit deterministic expectations
- Expected values are derived analytically and are not computed by v1 helpers
- Perturbing an expected edge or split fraction makes the test fail clearly
- Any v0 comparison is labeled as regression coverage for [kb/v0-errata/clock-min-dev-fixed-slope-score.md](../v0-errata/clock-min-dev-fixed-slope-score.md), not as a correctness oracle

## Locations

- Generic reroot implementation [packages/treetime/src/reroot/mod.rs](../../packages/treetime/src/reroot/mod.rs)
- Existing algebra cross-check [packages/treetime/src/reroot/__tests__/test_div_stats.rs](../../packages/treetime/src/reroot/__tests__/test_div_stats.rs)
- Clock dispatch coverage [packages/treetime/src/clock/__tests__/test_reroot.rs](../../packages/treetime/src/clock/__tests__/test_reroot.rs)
- v0 scoring defect [kb/v0-errata/clock-min-dev-fixed-slope-score.md](../v0-errata/clock-min-dev-fixed-slope-score.md)

## Related issues

- [N-reroot-split-optimizer-default-diverges-from-v0.md](N-reroot-split-optimizer-default-diverges-from-v0.md)
