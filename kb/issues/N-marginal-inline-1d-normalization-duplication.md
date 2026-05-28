# Inline 1D normalization duplicated in sparse marginal passes

Two adjacent code blocks in `marginal_passes.rs` perform identical 6-line sum-and-divide normalization with the same `norm > 0.0 && norm.is_finite()` guard and uniform-profile fallback. This is the 1D specialization of the 2D `normalize_inplace` in `marginal_core.rs`.

## Locations

- `packages/treetime/src/partition/marginal_passes.rs:402-409:` variable-sites normalization loop
- `packages/treetime/src/partition/marginal_passes.rs:416-423:` fixed-sites normalization loop

Both blocks:

```rust
let norm = dis.sum();
if norm > 0.0 && norm.is_finite() {
  delta_ll += norm.ln();
  dis / norm
} else {
  delta_ll += NEG_INFINITY;
  uniform
}
```

## Impact

If the normalization logic changes (e.g., handling of subnormals, log-space conversion), both blocks need identical updates.

## Action

Extract a `normalize_1d_inplace(dis, uniform) -> (Array1<f64>, f64)` helper that returns the normalized distribution and log-likelihood contribution.
