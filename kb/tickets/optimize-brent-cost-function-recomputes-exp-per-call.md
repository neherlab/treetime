# Brent cost function recomputes exp(eigvals\*t) per call

When multiple `OptimizationContribution` entries share the same GTR eigvals at the same branch length, the `exp_ev` vector is redundantly computed per contribution.

v1: [`packages/treetime/src/optimize/method_brent.rs`](../../packages/treetime/src/optimize/method_brent.rs)

## Related issues

- Source: [kb/issues/N-code-quality-conventions.md](../issues/N-code-quality-conventions.md) -- delete after full resolution
