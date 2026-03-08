# Dense backward pass produces NaN for all-zero probability rows

`normalize_from_log()` (`#normalize_from_log`) in
[`packages/treetime/src/representation/partition/marginal_dense.rs#L318-L345`](../../packages/treetime/src/representation/partition/marginal_dense.rs#L318-L345)
produces NaN when all log-probabilities in a row are `-inf` (all states have
zero probability at a position).

## Mechanism

The dense backward pass combines child messages in log space at
[`marginal_dense.rs#L167`](../../packages/treetime/src/representation/partition/marginal_dense.rs#L167):

```rust
let mut log_dis = msgs[0].mapv(f64::ln);
```

If any probability is exactly `0.0`, `f64::ln(0.0) = -inf`. When ALL states at
a position have zero probability (all entries `-inf`), `normalize_from_log`
computes:

```
max_val = -inf
(log_val - max_val).exp() = (-inf - (-inf)).exp() = NaN
```

IEEE 754: `-inf - (-inf) = NaN`. The NaN propagates through row normalization
and log-likelihood accumulation, silently corrupting downstream results.

## Trigger conditions

All states at a position must have zero probability after combining child
messages. This requires every child to assign zero probability to every state
at that position. Plausible scenarios:

- Degenerate GTR models with zero off-diagonal entries
- Extreme branch lengths causing numerical underflow in `expQt` before the
  log conversion
- Corrupted or pathological input alignments

In typical operation with well-conditioned GTR models and moderate branch
lengths, at least one state retains positive probability at every position.

## Fix direction

Guard the all-`-inf` case in `normalize_from_log`: when `max_val` is
`-inf`, set the row to uniform distribution (or zero with log-likelihood
penalty) instead of computing `NaN`.
