# Skyline grid construction accepts invalid endpoint arrays

> **Partially obsolete.** The `SkylineCostFunction` and the externally supplied
> `log_tc` array are gone: the segment grid is now derived internally as equal-width
> boundaries over the tree's time span (`equal_width_boundaries()`), and the
> optimizer is a convex Newton solve. The remaining validation gap is narrower —
> internal `boundaries` indexing still lacks a typed nonempty/ordering guard. See
> [decisions/coalescent-skyline-convex-log-tc.md](../decisions/coalescent-skyline-convex-log-tc.md).
> The original text below predates the rewrite.

Skyline construction indexes the first and last time-grid values before a typed boundary proves that the time and log-$T_c$ arrays are nonempty, equal in length, finite, and strictly ordered. Objective tests also exercise reconstructed formula fragments instead of the production cost function.

## Potential solutions

- O1. Parse both arrays into one validated skyline-grid type whose constructor establishes every length, finiteness, and ordering invariant.
- O2. Repeat guards in each constructor and objective entry point. This can reject invalid arrays but permits the guards to drift and leaves invalid intermediate states representable.

## Recommendation

Use O1: parse both arrays into a validated skyline-grid type before endpoint indexing, and test the production objective directly. This validation is independent of the open extrapolation, quadrature, and simplex policies.

## Related issues

- [N-coalescent-skyline-extrapolation-policy-undecided.md](N-coalescent-skyline-extrapolation-policy-undecided.md)
