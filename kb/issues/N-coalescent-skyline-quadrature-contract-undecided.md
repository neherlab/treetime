# Coalescent quadrature over the expected lineage counts lacks an accuracy contract

> [!IMPORTANT]
> **Decision required.** The integration method and its error contract for the integrals over the expected lineage counts are undecided, and the method needs an explicit parity disposition against the different v0 integration rule. The options are under "Design axis: integration method".

The coalescent statistics are integrals over the expected lineage counts: the pair-rate integral
$I_i = \int_i \tfrac12(\bar k^2 - \bar k + V)\,dt$ and, through the prior, $\int \bar k(t)/T_c(t)\,dt$
(`LineageExpectation::segment_statistics` [packages/treetime/src/coalescent/mean_field.rs#L374](../../packages/treetime/src/coalescent/mean_field.rs#L374)
and `CoalescentPrior::new` [packages/treetime/src/coalescent/mean_field.rs#L38](../../packages/treetime/src/coalescent/mean_field.rs#L38)). Both use the
trapezoid rule on a uniform grid of 4001 points (`EXPECTATION_GRID_POINTS` [packages/treetime/src/coalescent/mean_field.rs#L19](../../packages/treetime/src/coalescent/mean_field.rs#L19)) over the span of the node distributions, with $T_c$ evaluated at cell
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
