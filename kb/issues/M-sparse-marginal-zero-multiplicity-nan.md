# Sparse marginal zero-multiplicity fixed weight produces NaN

## Summary

In `fn combine_messages()`, multiplying a zero fixed-state weight by a `NEG_INFINITY` log-normalization produces NaN, corrupting the `msg_to_child.log_lh` and downstream branch-length distributions.

## Details

`fn combine_messages()` [`packages/treetime/src/partition/marginal/sparse/message.rs#L21-L126`](../../packages/treetime/src/partition/marginal/sparse/message.rs#L21-L126):

The line `seq_dis.log_lh += LogLh::new(weights[index] * log_norm)` [`packages/treetime/src/partition/marginal/sparse/message.rs#L119`](../../packages/treetime/src/partition/marginal/sparse/message.rs#L119) is unguarded. `weights[index]` is the composition count of a canonical state minus the candidate sites of that state. When it is 0 and `log_norm` is `f64::NEG_INFINITY` (from a zero-sum profile), the result is `0.0 * NEG_INFINITY = NaN`, although a state with zero multiplicity contributes exactly zero to the summed log likelihood. This NaN propagates through `msg_to_child.log_lh` into branch-length optimization.

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

## Proposed implementation plan

The cited line is current: [packages/treetime/src/partition/marginal/sparse/message.rs#L119](../../packages/treetime/src/partition/marginal/sparse/message.rs#L119) still adds `weights[index] * log_norm` without a guard; `weights[index]` is decremented per candidate site at `message.rs#L99`.

1. In `fn combine_messages()`, skip the term when `weights[index] == 0.0` (for example `filter` on the `izip!` iterator at `message.rs#L118`), and keep the product for positive weights
2. Optionally add a debug assertion that `weights[index] >= 0.0`, because a negative weight means more candidate sites than the composition count

Tests (`packages/treetime/src/partition/marginal/sparse/__tests__/`):

- Unit test of `combine_messages()` with a canonical state of zero multiplicity whose shared profile sums to zero: `log_lh` is finite and equals the sum of the other states' terms. Oracle: hand-computed sum (a zero-multiplicity term is the additive identity)
- Regression test on a rerooted sparse tree that leaves a canonical state with zero count: the branch log-likelihoods are finite and equal to dense within `1e-10`. Oracle: dense inference

Validation: `just test-rs`, `just lint-rs`.
