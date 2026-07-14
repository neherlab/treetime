# Sparse marginal fixed_counts zero-multiplicity produces NaN

## Summary

In `fn combine_messages()`, multiplying zero `fixed_counts` by `NEG_INFINITY` log-normalization produces NaN, corrupting the `msg_to_child.log_lh` and downstream branch-length distributions.

## Details

`packages/treetime/src/partition/marginal_helpers.rs`:

The line `seq_dis.log_lh += fixed_counts[&state] * log_norm` is unguarded. When `fixed_counts` is 0 for a character state and `log_norm` is `f64::NEG_INFINITY` (from a zero-sum profile), the result is `0.0 * NEG_INFINITY = NaN`. This NaN propagates through `msg_to_child.log_lh` into branch-length optimization.

Zero multiplicities are ordinary because canonical composition states are initialized even when their counts are zero. The defective combination requires a zero multiplicity paired with an all-negative-infinity normalization.

## Required behavior

Skip the likelihood term only when its multiplicity is zero. Preserve existing positive-multiplicity behavior; the separate normalization issue tracks its error contract.

## v0 comparison

v0 uses dense arrays and numpy operations that naturally handle zero multiplication without NaN (numpy's `0 * -inf = nan` but v0's flow avoids the zero-count case by construction since composition is always derived from a full sequence).

## Related tickets

- [kb/tickets/sparse-guard-zero-multiplicity-log-normalization.md](../tickets/sparse-guard-zero-multiplicity-log-normalization.md)
