# Dense-sparse duplication in OptimizationContribution enum dispatch

The eval loop duplication (4 identical inner loops across dense/sparse eval files) is resolved: `evaluate_site_contributions()` in `optimize_eval.rs` provides a single implementation, with thin wrappers in `optimize_dense_eval.rs` and `optimize_sparse_eval.rs`.

The same Dense-vs-Sparse iteration split still exists in `OptimizationContribution` methods in [packages/treetime/src/commands/optimize/optimize_unified.rs#L112-L184](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L112-L184):

| Method                             | Arms identical? | Notes                                                                  |
| :--------------------------------- | :-------------: | :--------------------------------------------------------------------- |
| `all_sites_valid_at_zero()`        |      near       | same `site_lh > 0.0 && site_lh.is_finite()` check, different iteration |
| `has_unimodal_branch_likelihood()` |     **yes**     | both arms access `.gtr.unimodal_branch_likelihood`                     |
| `zero_branch_length_derivative()`  |      near       | same weighted-eigenvalue formula, different iteration                  |

## Fix

Add a shared site-iteration accessor to `OptimizationContribution` returning `impl Iterator<Item = (f64, ArrayView1<f64>)>` (same interface as `evaluate_site_contributions`). The three methods above can then use the shared iterator, and `has_unimodal_branch_likelihood()` can access `.gtr` via a shared accessor.

## Impact

Low. Maintenance cost only, no correctness issue.
