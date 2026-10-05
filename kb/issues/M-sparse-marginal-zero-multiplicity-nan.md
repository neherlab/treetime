# Sparse marginal fixed_counts zero-multiplicity produces NaN

## Summary

In `fn combine_messages()`, multiplying zero `fixed_counts` by `NEG_INFINITY` log-normalization produces NaN, corrupting the `msg_to_child.log_lh` and downstream branch-length distributions.

## Details

`fn combine_messages()` [`packages/treetime/src/partition/marginal/sparse/message.rs#L22-L104`](../../packages/treetime/src/partition/marginal/sparse/message.rs#L22-L104):

The line `seq_dis.log_lh += LogLh::new(fixed_counts[&state] * log_norm)` [`packages/treetime/src/partition/marginal/sparse/message.rs#L99`](../../packages/treetime/src/partition/marginal/sparse/message.rs#L99) is unguarded. When `fixed_counts` is 0 for a character state and `log_norm` is `f64::NEG_INFINITY` (from a zero-sum profile), the result is `0.0 * NEG_INFINITY = NaN`, although a state with zero multiplicity contributes exactly zero to the summed log likelihood. This NaN propagates through `msg_to_child.log_lh` into branch-length optimization.

Zero multiplicities are ordinary because canonical composition states are initialized even when their counts are zero. The defective combination requires a zero multiplicity paired with an all-negative-infinity normalization.

## Required behavior

Skip the likelihood term only when its multiplicity is zero. Preserve existing positive-multiplicity behavior; the separate normalization issue tracks its error contract, and this fix does not broaden into it.

## Validation

- A regression test with zero multiplicity and negative-infinite normalization asserts a finite aggregate, without changing the contributions of other states
- A test exercises the rerooted sparse path that can produce the zero-count composition
- Zero-multiplicity terms act as additive identities and never introduce NaN
- Positive-multiplicity behavior is unchanged; sparse marginal and reroot regression tests pass

## v0 comparison

v0 uses dense arrays and numpy operations that naturally handle zero multiplication without NaN (numpy's `0 * -inf = nan` but v0's flow avoids the zero-count case by construction since composition is always derived from a full sequence).
