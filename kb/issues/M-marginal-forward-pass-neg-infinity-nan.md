# Marginal forward pass NEG_INFINITY produces NaN

## Summary

The marginal forward pass computes NaN from infinity arithmetic and poisons likelihood accumulators when nodes have degenerate (zero-probability) profiles.

## Details

Two related defects in the forward pass:

### NaN from log_lh subtraction

`packages/treetime/src/partition/marginal_passes.rs:355:` and `marginal_dense.rs:327:`

The forward pass computes:

```rust
log_lh: seq_info.profile.log_lh - child_edge_data.msg_from_child.log_lh
```

When both values are `f64::NEG_INFINITY` (degenerate node), the subtraction produces NaN (`-inf - (-inf) = NaN`). No guard prevents this arithmetic.

### NEG_INFINITY poisons delta_ll accumulator

`packages/treetime/src/partition/marginal_passes.rs:406:`

The sparse forward pass adds `NEG_INFINITY` directly into the `delta_ll` accumulator. A single degenerate site contribution poisons the `log_lh` for the entire child message, making the total log-likelihood meaningless for downstream convergence checks.

## Impact

- NaN propagates silently through the tree during forward pass
- Convergence metrics become unreliable
- Affects trees with degenerate nodes (zero-probability states at any position)

## Fix

Guard the subtraction: when both operands are `NEG_INFINITY`, treat the result as 0.0 (the normalization constants cancel). For the accumulator, skip `NEG_INFINITY` contributions or use a sentinel that does not poison finite sums.
