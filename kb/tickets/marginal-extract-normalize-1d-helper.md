# Extract 1D normalization helper from sparse marginal passes

Two adjacent code blocks in marginal_passes.rs perform identical 6-line sum-and-divide normalization.

## Current state

`partition/marginal_passes.rs:402-409` (variable sites) and `:416-423` (fixed sites) each inline the same normalization: compute sum, check `> 0.0 && is_finite()`, divide or fall back to uniform, accumulate log-likelihood.

## Target state

A `normalize_1d(dis: ArrayView1<f64>, uniform: &Array1<f64>) -> (Array1<f64>, f64)` helper returns the normalized distribution and the log-likelihood contribution (`norm.ln()` or `NEG_INFINITY`).

## Implementation

1. Define `normalize_1d` in `partition/marginal_passes.rs` or `partition/marginal_core.rs` near the existing 2D `normalize_inplace`
2. Replace both inline blocks with calls to the helper
3. Verify sparse marginal pass tests produce identical results

## Related issues

Source: `kb/issues/N-marginal-inline-1d-normalization-duplication.md`
