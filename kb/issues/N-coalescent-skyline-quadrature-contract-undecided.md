# Coalescent quadrature over the expected lineage counts lacks an accuracy contract

> [!IMPORTANT]
> **Decision required.** The integration method and its error contract for the integrals over the expected lineage counts are undecided, and the method needs an explicit parity disposition against the different v0 integration rule. The options are under "Design axis: integration method".

The coalescent statistics are integrals over the expected lineage counts: the pair-rate integral
$I_i = \int_i \tfrac12(\bar k^2 - \bar k + V)\,dt$ and, through the prior, $\int \bar k(t)/T_c(t)\,dt$
(`LineageExpectation::segment_statistics` and `CoalescentPrior::new` in
[packages/treetime/src/coalescent/mean_field.rs](../../packages/treetime/src/coalescent/mean_field.rs)). Both use the
trapezoid rule on a uniform grid of 4001 points over the span of the node distributions, with $T_c$ evaluated at cell
midpoints. The rule has no stated error bound: a step of $\bar k$ at point node times costs half a grid cell (the
tests allow a relative error of $10^{-4}$ on $T_c = I/M$ for that reason), and the grid spacing scales with the span of
the deepest node distribution rather than with the variation of $\bar k$ or $T_c$. Since phase 2 of
[timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md)
the removed midpoint quadrature of the point-time lineage counts (`compute_integral_merger_rate`) no longer applies.
V0 uses a different numerical integration rule, so the selected algorithm also requires an explicit parity
disposition.

The numerical method must define an observable error contract for these integrals over each segment.

## Design axis: integration method

- O1. Retain midpoint integration with adaptive interval refinement and a declared absolute/relative error criterion.
- O2. Use an adaptive maintained quadrature implementation directly over the piecewise skyline.
- O3. Derive an exact interval integral for the selected interpolation of $T_c(t)$, with stable handling of limiting cases.

## Recommendation

Prefer O3 when the selected interpolation has a stable analytic integral; otherwise use O2. Validate convergence against an independent high-precision integral before implementation starts.

## Related issues
