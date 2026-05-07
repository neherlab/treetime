# Promote debug_assert to production guards for gamma and branch likelihood

`debug_assert!(gamma > 0.0)` is stripped in release builds. When gamma <= 0 reaches this code path in production, the computation produces inf/NaN silently.

v1: [`packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L61`](../../packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L61)

## Impact

Silent inf/NaN propagation in production builds for edge-case inputs.

## Related issues

- Source: [N-numerical-stability-magic-constants.md](../issues/N-numerical-stability-magic-constants.md) -- delete after full resolution
- See also: [N-error-suppression-unwrap-or-defaults.md](../issues/N-error-suppression-unwrap-or-defaults.md) (related debug_assert patterns)
