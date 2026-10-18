# Min-dev rerooting lacks oracles for endpoint roots and ties

## Background

The intended `min-dev` objective minimizes root-to-tip distance variance. Clock `RerootMethod::MinDev` and `optimize --reroot=min-dev` run the same clock root search at `RootObjective::FixedRate(0.0)`, exact on every edge ([kb/algo/reroot.md](../algo/reroot.md)). Analytic fixtures cover interior optima: a cherry roots at equal distances and the star $(A{:}1, B{:}2, C{:}6)$ roots $9/4$ up the branch to $C$ ([packages/treetime/src/reroot/__tests__/test_orchestrate.rs](../../packages/treetime/src/reroot/__tests__/test_orchestrate.rs)); the clock dispatch is independent of the date values ([packages/treetime/src/clock/__tests__/test_reroot.rs](../../packages/treetime/src/clock/__tests__/test_reroot.rs)).

TreeTime v0 cannot serve as the correctness oracle for this objective. Although v0 passes `slope=0` and documents minimum-variance behavior, its candidate score retains the estimated-rate regression term. See [kb/v0-errata/clock-min-dev-fixed-slope-score.md](../v0-errata/clock-min-dev-fixed-slope-score.md).

## Impact

Negligible. Interior optima are tested analytically; the residual risk lies in endpoint handling (an optimum at an existing node, clamped splits) and in ties between edges.

## Required oracle

- Fixtures whose minimum-variance root lies on an existing node, with the expected node derived analytically
- Fixtures with tied optima on several edges, with an explicit deterministic expectation
- Any v0 comparison is labeled as regression coverage for [kb/v0-errata/clock-min-dev-fixed-slope-score.md](../v0-errata/clock-min-dev-fixed-slope-score.md), not as a correctness oracle
