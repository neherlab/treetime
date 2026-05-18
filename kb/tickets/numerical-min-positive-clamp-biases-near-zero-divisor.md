# Discrete MIN_POSITIVE clamp biases near-zero divisor

Silently biases near-zero values to `f64::MIN_POSITIVE` without propagating information about the degenerate log-likelihood to callers.

v1: [`packages/treetime/src/partition/discrete.rs#L207`](../../packages/treetime/src/partition/discrete.rs#L207)

## Impact

Degenerate inputs silently clamped rather than rejected or surfaced to callers.

## Related issues

- Source: [N-numerical-stability-magic-constants.md](../issues/N-numerical-stability-magic-constants.md) -- delete after full resolution
