# Shared marginal forward pass MIN_POSITIVE clamp biases near-zero divisor

Silently biases near-zero values to `f64::MIN_POSITIVE` without propagating information about the degenerate log-likelihood to callers. Used by both dense and discrete partitions.

v1: [`packages/treetime/src/partition/marginal_core.rs#L163`](../../packages/treetime/src/partition/marginal_core.rs#L163)

## Impact

Degenerate inputs silently clamped rather than rejected or surfaced to callers.

## Related issues

- Source: [N-numerical-stability-magic-constants.md](../issues/N-numerical-stability-magic-constants.md)
- Source: [N-marginal-forward-zero-divisor-floor.md](../issues/N-marginal-forward-zero-divisor-floor.md)
