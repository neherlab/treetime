# Unguarded ln(site_lh) in evaluation functions

The evaluation functions `evaluate_dense_contribution_impl()` and `evaluate_sparse_contribution_impl()` compute `site_lh.ln()` without checking `site_lh > 0`. When per-site coefficients sum to zero or negative (degenerate eigenvalue-space projections), `ln(0)` produces `-inf` and `ln(negative)` produces `NaN`. These propagate silently through the Newton step and can corrupt branch length estimates.

## Location

- [packages/treetime/src/commands/optimize/optimize_dense_eval.rs#L33](../../packages/treetime/src/commands/optimize/optimize_dense_eval.rs#L33): `log_lh += site_lh.ln();`
- [packages/treetime/src/commands/optimize/optimize_dense_eval.rs#L41](../../packages/treetime/src/commands/optimize/optimize_dense_eval.rs#L41): same in no-derivatives path
- [packages/treetime/src/commands/optimize/optimize_sparse_eval.rs#L34](../../packages/treetime/src/commands/optimize/optimize_sparse_eval.rs#L34): `log_lh += multiplicity * site_lh.ln();`
- [packages/treetime/src/commands/optimize/optimize_sparse_eval.rs#L46](../../packages/treetime/src/commands/optimize/optimize_sparse_eval.rs#L46): same in no-derivatives path

## Partial mitigation

`is_zero_branch_optimal()` checks `all_sites_valid_at_zero()` which verifies `site_lh > 0 && site_lh.is_finite()` before the zero-branch shortcut. But the main Newton and grid search evaluation paths do not guard against this. For valid phylogenetic data, coefficients produce positive site likelihoods. The risk is degenerate input or numerical edge cases.

## Impact

Medium. NaN/Inf corruption would produce silently wrong branch lengths. Mitigated by the fact that valid phylogenetic inputs produce positive coefficient sums.

## Fix

Add `debug_assert!(site_lh > 0.0)` or return early with a sentinel value when `site_lh <= 0.0`.
