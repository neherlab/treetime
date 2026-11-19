# Coalescent quadrature over the expected lineage counts lacks an accuracy contract

> [!IMPORTANT]
> **Decision required.** The integration method and its error contract for the integrals over the expected lineage counts are undecided, and the method needs an explicit parity disposition against the different v0 integration rule. The options are under "Design axis: integration method".

The coalescent statistics are integrals over the expected lineage counts: the pair-rate integral
$I_i = \int_i \tfrac12(\bar k^2 - \bar k + V)\,dt$ and, through the prior, $\int \bar k(t)/T_c(t)\,dt$
(`LineageExpectation::segment_statistics` [packages/treetime/src/coalescent/mean_field.rs#L484](../../packages/treetime/src/coalescent/mean_field.rs#L484)
and `CoalescentPrior::new` [packages/treetime/src/coalescent/mean_field.rs#L42](../../packages/treetime/src/coalescent/mean_field.rs#L42)). Both use the
trapezoid rule on a uniform grid of 4001 points (`EXPECTATION_GRID_POINTS` [packages/treetime/src/coalescent/mean_field.rs#L22](../../packages/treetime/src/coalescent/mean_field.rs#L22)) over the span of the node distributions, with $T_c$ evaluated at cell
midpoints. The steps of $\bar k$ at point node times are exact: each step time enters the grid as a zero-width pair
of points that holds the left and the right limit, and cells are split at segment boundaries. The rest of the rule has
no stated error bound: the grid spacing scales with the span of the deepest node distribution rather than with the
variation of the smooth part of $\bar k$ or of $T_c$, and $T_c$ is evaluated at cell midpoints.
V0 uses a different numerical integration rule, so the selected algorithm also requires an explicit parity
disposition.

The numerical method must define an observable error contract for these integrals over each segment.

## Design axis: integration method

- O1. Retain the current rule (trapezoid rule in $\bar k$, $T_c$ at cell midpoints) with adaptive interval refinement and a declared absolute/relative error criterion.
- O2. Use an adaptive maintained quadrature implementation directly over the piecewise skyline.
- O3. Derive an exact interval integral for the selected interpolation of $T_c(t)$, with stable handling of limiting cases.

## Recommendation

Prefer O3 when the selected interpolation has a stable analytic integral; otherwise use O2. Validate convergence against an independent high-precision integral before implementation starts.

## Proposed implementation plan

The locations are current (`EXPECTATION_GRID_POINTS` at `mean_field.rs` L22, `CoalescentPrior::new` at L42, `LineageExpectation::segment_statistics` at L494). The plan depends on the open decision.

1. Before the choice: measure the error of the current rule against an independent high-precision integral (`mpmath` quadrature of the same $\bar k(t)$, $V(t)$ and $T_c(t)$, exported from a small tree) for each segment, on a tree with point dates and on one with wide date ranges; report the error per segment and its scaling with the grid size
2. **O1**: refine grid cells where the difference between the trapezoid and the midpoint estimates of a cell exceeds the declared tolerance, in `segment_statistics` and `CoalescentPrior::new`
3. **O2**: replace the trapezoid sums with a maintained adaptive quadrature (vetted through `dependencies.instructions.md`) over each skyline segment, with the step times of $\bar k$ as breakpoints
4. **O3**: integrate each cell exactly for the piecewise-constant $T_c$ of the skyline (`kb/decisions/coalescent-skyline-convex-log-tc.md`) and the selected interpolation of $\bar k$
5. Record the chosen method and its parity disposition against v0 in `kb/decisions/`
6. Tests in `packages/treetime/src/coalescent/__tests__/`: per-segment integrals on a small tree within the declared error criterion (oracle: the `mpmath` reference values of step 1, cited at the expected-value declaration); a pure step function $\bar k$ with constant $T_c$ integrates exactly (oracle: manual calculation)
7. Validation: `just check`, `just test-rs`
