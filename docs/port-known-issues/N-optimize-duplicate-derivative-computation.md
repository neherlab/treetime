# Derivative ratio computed twice per site in evaluation functions

The per-site first derivative ratio `(coefficients * &ev_exp_ev).sum() / site_lh` is computed twice in the Newton evaluation path: once for the first derivative and once for the squared term in the second derivative. Caching the ratio eliminates one temporary allocation and element-wise multiply per site per edge per Newton iteration.

## Location

- [packages/treetime/src/commands/optimize/optimize_dense_eval.rs#L34-L36](../../packages/treetime/src/commands/optimize/optimize_dense_eval.rs#L34-L36)
- [packages/treetime/src/commands/optimize/optimize_sparse_eval.rs#L35-L37](../../packages/treetime/src/commands/optimize/optimize_sparse_eval.rs#L35-L37)
- [packages/treetime/src/commands/optimize/optimize_dense.rs#L52-L53](../../packages/treetime/src/commands/optimize/optimize_dense.rs#L52-L53)

```rust
derivative += (coefficients * &ev_exp_ev).sum() / site_lh;
second_derivative += (coefficients * &ev2_exp_ev).sum() / site_lh
  - ((coefficients * &ev_exp_ev).sum() / site_lh).powi(2);  // recomputed
```

## Impact

Negligible for small alignments. For large alignments with many sites per edge, the duplicate computation doubles the per-site work in the derivative path. The allocation is a 4-element (nucleotide) or 20-element (amino acid) temporary array.

## Fix

Cache the ratio:

```rust
let d1 = (coefficients * &ev_exp_ev).sum() / site_lh;
derivative += d1;
second_derivative += (coefficients * &ev2_exp_ev).sum() / site_lh - d1.powi(2);
```

This also makes the Hessian formula clearer and prevents the multiplicity-squaring bug pattern.
