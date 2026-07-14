# Prevent zero-multiplicity sparse likelihood terms from producing NaN

`fn combine_messages()` [`packages/treetime/src/partition/marginal_helpers.rs#L30-L84`](../../packages/treetime/src/partition/marginal_helpers.rs#L30-L84) multiplies a fixed-state multiplicity by a log-normalization term. In Rust, `0.0 * f64::NEG_INFINITY` evaluates to NaN even though a state with zero multiplicity contributes exactly zero to the summed log likelihood.

## Implementation

1. In `fn combine_messages()` [`packages/treetime/src/partition/marginal_helpers.rs#L30-L84`](../../packages/treetime/src/partition/marginal_helpers.rs#L30-L84), skip the likelihood term when its multiplicity is zero.
2. Leave positive-multiplicity behavior unchanged; do not broaden this fix into the separate normalization contract.
3. Add a regression with zero multiplicity and negative-infinite normalization that asserts a finite aggregate without changing other state contributions.
4. Exercise the rerooted sparse path that can produce the reported zero-count composition.

## Acceptance criteria

- Zero-multiplicity terms are additive identities and never introduce NaN.
- Positive-multiplicity behavior is unchanged.
- Sparse marginal and reroot regression tests pass.

## Related issues

- Source: [kb/issues/M-sparse-marginal-zero-multiplicity-nan.md](../issues/M-sparse-marginal-zero-multiplicity-nan.md)
