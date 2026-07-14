# Constant-$T_c$ optimization can report success without a parameter

When optimization returns no best parameter, the implementation falls back to the initial $T_c$ value. Success is then based on finiteness alone, so an optimization that failed to produce a parameter can be reported as successful.

`fn optimize_tc()` applies `unwrap_or(initial_log_tc)` to `best_param` and derives `success` only from finite parameter and cost values [packages/treetime/src/coalescent/optimize_tc.rs#L55-L90](../../packages/treetime/src/coalescent/optimize_tc.rs#L55-L90).

## Potential solutions

- O1. Make successful optimization contain a validated parameter and termination record.
- O2. Return the initial parameter with a distinct non-success status. This preserves a usable fallback without misreporting optimization.

## Recommendation

Require a successful optimizer termination and a returned, finite, domain-valid parameter. Preserve the optimizer’s original diagnostic on failure; do not replace absence with the initial value in a success result.

## Related issues

- [N-timetree-negative-coalescent-tc.md](N-timetree-negative-coalescent-tc.md)
- [N-coalescent-skyline-simplex-initialization-undecided.md](N-coalescent-skyline-simplex-initialization-undecided.md)
