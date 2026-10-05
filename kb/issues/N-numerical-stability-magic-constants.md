# Numerical stability magic constants

## Summary

Numeric thresholds remain hardcoded without named constants, documentation, or a defined fallback contract for degenerate inputs.

## Instances

### 1e-10 magic denominator in relaxed_clock.rs

`fn apply_relaxed_clock()` uses `1e-10` as a denominator floor in three places [packages/treetime/src/timetree/optimization/relaxed_clock.rs#L46](../../packages/treetime/src/timetree/optimization/relaxed_clock.rs#L46), [#L66](../../packages/treetime/src/timetree/optimization/relaxed_clock.rs#L66), [#L80](../../packages/treetime/src/timetree/optimization/relaxed_clock.rs#L80). The constant has no name. Below the floor, the code skips the child contribution or sets the rate multiplier `gamma` to `1.0` without a warning.

v0's `relaxed_clock()` has no such guard and divides directly [packages/legacy/treetime/treetime/treetime.py#L1110-L1130](../../packages/legacy/treetime/treetime/treetime.py#L1110-L1130), so the guard and its fallback values are an undocumented divergence from v0.

## Impact

- Degenerate relaxed-clock coefficients silently produce default rate multipliers

## Related issues

- [N-marginal-forward-zero-divisor-floor.md](N-marginal-forward-zero-divisor-floor.md): the `f64::MIN_POSITIVE` divisor floor of the marginal forward pass [packages/treetime/src/partition/marginal/shared/normalize.rs#L59](../../packages/treetime/src/partition/marginal/shared/normalize.rs#L59)
- [N-timetree-branch-likelihood-positivity-checked-only-in-debug.md](N-timetree-branch-likelihood-positivity-checked-only-in-debug.md)
