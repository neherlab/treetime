# expQt uses O(n^3) diagonal matrix multiply instead of O(n^2) broadcast

`GTR::expQt()` constructs a full diagonal matrix via `Array2::from_diag()` then multiplies with two `dot()` calls, producing O(n^3) work. A row-wise broadcast scaling achieves the same result in O(n^2).

## Location

- [packages/treetime/src/gtr/gtr.rs#L471-L480](../../packages/treetime/src/gtr/gtr.rs#L471-L480): `expQt()` and `exp_lt()`
- [packages/treetime/src/gtr/gtr.rs#L329-L333](../../packages/treetime/src/gtr/gtr.rs#L329-L333): `expQt_with_rate()`

```rust
let eLambdaT = Array2::from_diag(&self.exp_lt(t));
let Qt = self.v.dot(&eLambdaT.dot(&self.v_inv));
```

## Impact

Negligible for nucleotides (n=4, 64 vs 16 ops). Measurable for amino acids (n=20, 8000 vs 400 ops per call). `expQt` is called once per edge per marginal reconstruction pass.

## Fix

Replace `from_diag` + dot with broadcast row scaling:

```rust
let exp_ev = self.exp_lt(t);
let scaled_v_inv = &self.v_inv * &exp_ev.insert_axis(Axis(1));
let Qt = self.v.dot(&scaled_v_inv);
```

Or use `v * exp_ev` (column scaling) depending on the layout convention.
