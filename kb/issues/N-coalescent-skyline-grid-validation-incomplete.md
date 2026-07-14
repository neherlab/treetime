# Skyline grid construction accepts invalid endpoint arrays

Skyline construction indexes the first and last time-grid values before a typed boundary proves that the time and log-$T_c$ arrays are nonempty, equal in length, finite, and strictly ordered [packages/treetime/src/coalescent/skyline.rs#L215-L227](../../packages/treetime/src/coalescent/skyline.rs#L215-L227). Objective tests also exercise reconstructed formula fragments instead of the production `impl CostFunction for &SkylineCostFunction` [packages/treetime/src/coalescent/skyline.rs#L160-L200](../../packages/treetime/src/coalescent/skyline.rs#L160-L200).

## Potential solutions

- O1. Parse both arrays into one validated skyline-grid type whose constructor establishes every length, finiteness, and ordering invariant.
- O2. Repeat guards in each constructor and objective entry point. This can reject invalid arrays but permits the guards to drift and leaves invalid intermediate states representable.

## Recommendation

Use O1: parse both arrays into a validated skyline-grid type before endpoint indexing, and test the production objective directly. This validation is independent of the open extrapolation, quadrature, and simplex policies.

## Related issues

- [N-coalescent-skyline-extrapolation-policy-undecided.md](N-coalescent-skyline-extrapolation-policy-undecided.md)
