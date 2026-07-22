# Skyline merger-rate quadrature lacks an accuracy contract

> **Scope narrowed.** The skyline is now piecewise-*constant*, and its optimizer
> attributes each lineage-count interval wholly to the segment of its midpoint —
> the same midpoint convention the model uses — so the midpoint rule is *exact* for
> the piecewise-constant skyline (no quadrature error, and the analytic optimum
> matches the model-evaluated likelihood). See
> [decisions/coalescent-skyline-convex-inverse-tc.md](../decisions/coalescent-skyline-convex-inverse-tc.md).
> The contract below still applies to genuinely continuous $T_c(t)$ (e.g. a
> `Formula` distribution), which the ignored smooth-$T_c$ test exercises.

`fn compute_integral_merger_rate()` evaluates $T_c(t)$ once at each lineage-count interval midpoint [packages/treetime/src/coalescent/integration.rs#L44-L71](../../packages/treetime/src/coalescent/integration.rs#L44-L71). The approximation has no stated convergence order, refinement rule, or error bound for a varying skyline.

The ignored analytical test integrates $1/(0.01+0.004t)$ on $[0,10]$ and observes an error of approximately $1.6\times10^{-4}$ [packages/treetime/src/coalescent/__tests__/test_integration.rs#L157-L164](../../packages/treetime/src/coalescent/__tests__/test_integration.rs#L157-L164). V0 uses a different numerical integration rule, so the selected algorithm also requires an explicit parity disposition.

Let $\kappa(t)$ be the per-branch merger rate and let $I(t)$ be its cumulative integral:

$$I(t) = \int_0^t \kappa(s)\,ds.$$

The numerical method must define an observable error contract for approximating $I(t)$ over each interval.

## Design axis: integration method

- O1. Retain midpoint integration with adaptive interval refinement and a declared absolute/relative error criterion.
- O2. Use an adaptive maintained quadrature implementation directly over the piecewise skyline.
- O3. Derive an exact interval integral for the selected interpolation of $T_c(t)$, with stable handling of limiting cases.

## Recommendation

Prefer O3 when the selected interpolation has a stable analytic integral; otherwise use O2. Validate convergence against an independent high-precision integral before creating an implementation ticket.

## Related issues

- [N-coalescent-skyline-simplex-initialization-undecided.md](N-coalescent-skyline-simplex-initialization-undecided.md)
