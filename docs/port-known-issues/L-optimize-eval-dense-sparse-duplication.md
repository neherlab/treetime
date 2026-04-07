# Dense and sparse evaluation loops duplicated across 4 inner loops

The branch-length evaluation functions `evaluate_dense_contribution_impl()` and `evaluate_sparse_contribution_impl()` implement the same eigenvalue-space log-likelihood formula. Each function has an `if compute_derivatives` branch, giving 4 structurally identical inner loops across 2 files:

| File                                                                                                                                                     | Branch              | Multiplicity          |
| :------------------------------------------------------------------------------------------------------------------------------------------------------- | :------------------ | :-------------------- |
| [packages/treetime/src/commands/optimize/optimize_dense_eval.rs#L31-L37](../../packages/treetime/src/commands/optimize/optimize_dense_eval.rs#L31-L37)   | with derivatives    | implicit 1.0          |
| [packages/treetime/src/commands/optimize/optimize_dense_eval.rs#L39-L42](../../packages/treetime/src/commands/optimize/optimize_dense_eval.rs#L39-L42)   | without derivatives | implicit 1.0          |
| [packages/treetime/src/commands/optimize/optimize_sparse_eval.rs#L33-L38](../../packages/treetime/src/commands/optimize/optimize_sparse_eval.rs#L33-L38) | with derivatives    | explicit multiplicity |
| [packages/treetime/src/commands/optimize/optimize_sparse_eval.rs#L45-L47](../../packages/treetime/src/commands/optimize/optimize_sparse_eval.rs#L45-L47) | without derivatives | explicit multiplicity |

The only semantic difference between dense and sparse is iteration: dense iterates rows of `Array2<f64>` (each row is one site, implicit multiplicity 1.0), sparse iterates `Vec<SiteContribution>` (each element has explicit `multiplicity: f64` and `coefficients: Array1<f64>`). The per-site formula is identical.

## Evidence

Commit `ecea678d` ("fix(optimize): prevent ln(0) and division-by-zero at zero branch length") required changing the same `debug_assert!` in all 4 loops. Earlier commit `6765f079` ("perf(optimize): cache first-derivative ratio in evaluation loops") similarly required applying the same `d1` caching optimization in all 4 derivative loops. Any per-site change must be replicated 4 times.

## Impact

Low. Maintenance cost: formula changes, assertions, or instrumentation require 4 identical edits. No correctness issue, but the duplication is a defect magnet (the sparse Hessian multiplicity bug in `c56a0a8d` existed in one copy but not the other).

## Fix

Unify into a single evaluator taking `impl Iterator<Item = (f64, ArrayView1<f64>)>` (multiplicity, coefficients). Dense wraps `coefficients.outer_iter().map(|row| (1.0, row))`; sparse wraps `site_contributions.iter().map(|sc| (sc.multiplicity, sc.coefficients.view()))`. The 4 loops collapse to 2 (with/without derivatives), each written once.

Keep thin wrapper functions (`evaluate_dense_contribution`, `evaluate_sparse_contribution`) to preserve existing call sites.

Detailed implementation plan in [.memory/99-optimize-eval-dedup/refactor1.md](../../.memory/99-optimize-eval-dedup/refactor1.md).

## Downstream: `OptimizationContribution` enum dispatch

The same Dense-vs-Sparse iteration split propagates into `OptimizationContribution` methods in [packages/treetime/src/commands/optimize/optimize_unified.rs#L112-L184](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L112-L184):

| Method                             | Arms identical? | Notes                                                                       |
| :--------------------------------- | :-------------: | :-------------------------------------------------------------------------- |
| `evaluate()`                       |       no        | delegates to `evaluate_dense_contribution` / `evaluate_sparse_contribution` |
| `all_sites_valid_at_zero()`        |      near       | same `site_lh > 0.0 && site_lh.is_finite()` check, different iteration      |
| `has_unimodal_branch_likelihood()` |     **yes**     | both arms access `.gtr.unimodal_branch_likelihood`                          |
| `zero_branch_length_derivative()`  |      near       | same weighted-eigenvalue formula, different iteration                       |

Once the eval functions are unified behind a shared site iterator (see Fix above), `evaluate()` dispatch collapses to a single call. `all_sites_valid_at_zero()` and `zero_branch_length_derivative()` can use the same iterator. `has_unimodal_branch_likelihood()` can access `.gtr` directly via a shared accessor, eliminating the match entirely.

This is not a separate issue - it shares the same root cause (Dense and Sparse lack a common site-iteration interface) and the same fix resolves both.

## Related

- [Dead `optimize_dense::evaluate()` function](N-optimize-dense-evaluate-dead-code.md) - a third copy of the same formula, already dead code
