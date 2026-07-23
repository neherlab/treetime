# Skyline grid construction accepts invalid endpoint arrays

> **Validation added; only the typed-grid refactor remains.** `optimize_skyline`
> now rejects an over-segmented grid up front (requested segments exceeding the
> number of distinct merger times) and validates strictly increasing boundaries and
> positive per-segment exposure as backstops, erroring with an actionable message
> instead of masking the degeneracy with the pooled warm-start fallback
> ([packages/treetime/src/coalescent/skyline.rs](../../packages/treetime/src/coalescent/skyline.rs)).
> The skyline tests exercise the production `optimize_skyline` objective directly and
> assert the invariant and the over-segmentation error. The core defect (silent
> acceptance of invalid segmentations) is resolved; what remains is the optional O1
> refactor that moves these inline guards into a single validated skyline-grid type.
> See [decisions/coalescent-skyline-convex-inverse-tc.md](../decisions/coalescent-skyline-convex-inverse-tc.md).
> The original text below predates the rewrite.

Skyline construction indexes the first and last time-grid values before a typed boundary proves that the time and log-$T_c$ arrays are nonempty, equal in length, finite, and strictly ordered. Objective tests also exercise reconstructed formula fragments instead of the production cost function.

## Remaining work

- O1. Parse the segment grid into one validated skyline-grid type whose constructor establishes every length, finiteness, and ordering invariant, replacing the current inline guards. This is a structural cleanup; the invariants themselves are already enforced.

## Related issues

- [N-coalescent-skyline-extrapolation-policy-undecided.md](N-coalescent-skyline-extrapolation-policy-undecided.md)
