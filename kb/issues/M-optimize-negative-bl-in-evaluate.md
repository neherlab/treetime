# evaluate_site_contributions does not guard against negative branch lengths

`evaluate_site_contributions` at [packages/treetime/src/optimize/eval.rs#L42](../../packages/treetime/src/optimize/eval.rs#L42) computes `(eigvals * branch_length).mapv(f64::exp)`. For negative `branch_length` with negative eigenvalues, this computes `exp(positive * |t|)`, which grows exponentially and produces physically meaningless but finite likelihood values. The subsequent `.ln()` produces a finite log-likelihood that the optimizer treats as valid.

The `debug_assert!(site_lh.is_finite())` guard fires only in debug builds and does not catch all cases (the exponential overflow must exceed `f64::MAX` to produce `inf`).

## Impact

Timetree output with negative branch lengths (e.g., `bl = -1.89e-4`) silently flows into `exp(eigvals * t)` and then into the Brent cost function, producing meaningless optimization results without any diagnostic.

## Affected code

- Evaluation: [packages/treetime/src/optimize/eval.rs#L42](../../packages/treetime/src/optimize/eval.rs#L42)
- Related validation gap: [M-optimize-negative-branch-length-validation.md](M-optimize-negative-branch-length-validation.md)

## Fix

Add an explicit `branch_length >= 0` check at the entry to `evaluate_site_contributions`, returning an error for negative values. Callers already handle `Result` return types.

## Related

- [M-optimize-negative-branch-length-validation.md](M-optimize-negative-branch-length-validation.md): missing negative branch-length validation across optimizer modes
