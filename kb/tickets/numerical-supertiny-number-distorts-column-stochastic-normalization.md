# SUPERTINY_NUMBER distorts column-stochastic normalization

`SUPERTINY_NUMBER` (1e-24) added to `expQt` result distorts column-stochastic normalization. The additive constant prevents exact zeros but shifts column sums away from 1.0.

v1: [`packages/treetime/src/gtr/infer_gtr/dense.rs#L199`](../../packages/treetime/src/gtr/infer_gtr/dense.rs#L199)

## Impact

Column-stochastic property violated by additive perturbation.

## Related issues

- Source: [N-numerical-stability-magic-constants.md](../issues/N-numerical-stability-magic-constants.md) -- delete after full resolution
