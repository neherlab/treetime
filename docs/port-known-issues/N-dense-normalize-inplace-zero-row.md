# Dense normalize_inplace produces NaN for all-zero probability rows

`normalize_inplace()` (`#normalize_inplace`) in [`packages/treetime/src/representation/partition/marginal_dense.rs#L393-L399`](../../packages/treetime/src/representation/partition/marginal_dense.rs#L393-L399) divides each row by its sum. When a row sums to zero, `0.0 / 0.0 = NaN` (IEEE 754). The `0.0_f64.ln() = -inf` in the log-likelihood accumulation is correct, but the NaN in the output row propagates to downstream computation.

## Mechanism

```rust
fn normalize_inplace(dis: &mut Array2<f64>) -> f64 {
  let norm = dis.sum_axis(Axis(1));
  for (ri, mut row) in dis.outer_iter_mut().enumerate() {
    row /= norm[ri];  // NaN when norm[ri] == 0.0
  }
  norm.mapv(f64::ln).sum()
}
```

## Trigger conditions

All states at a position must have zero probability in the probability-space matrix. Plausible paths:

- Root node (line 278): `msg_to_parent.dis * gtr.pi` -- GTR equilibrium multiplication can underflow when both factors are near-zero
- Forward pass (lines 324, 335): products of parent and child distributions -- any factor with zeros produces all-zero products

The `normalize_from_log` all-`-inf` guard (which returns uniform instead of NaN) prevents one source of all-zero rows from reaching `normalize_inplace`. Other sources remain.

## Fix direction

Guard the zero-sum case: when `norm[ri] == 0.0`, set the row to uniform distribution, matching the `normalize_from_log` and sparse `logsumexp_normalize` convention.

## Related

- The stabilized log-space normalization path handles the same class of numerical degeneracy earlier in the pipeline. This issue covers the remaining probability-space path.
