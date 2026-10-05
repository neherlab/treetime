# Branch likelihood evaluation recomputes exp(eigvals\*t) per contribution

Every evaluation of the branch-length objective computes the vector $\exp(\lambda t)$ of the GTR eigenvalues $\lambda$ at branch length $t$ once for each `OptimizationContribution`. When several contributions on one edge share the same GTR eigenvalues, the same vector is computed several times per evaluation.

## Locations

- `fn evaluate_mixed_impl()` [packages/treetime/src/optimize/likelihood.rs#L44-L58](../../packages/treetime/src/optimize/likelihood.rs#L44-L58) loops over the contributions of an edge and evaluates each one separately
- `fn evaluate_dense_contribution()` and `fn evaluate_sparse_contribution()` pass the eigenvalues of the contribution's own GTR [packages/treetime/src/optimize/dense_eval.rs#L12](../../packages/treetime/src/optimize/dense_eval.rs#L12), [packages/treetime/src/optimize/sparse_eval.rs#L15](../../packages/treetime/src/optimize/sparse_eval.rs#L15)
- `fn evaluate_site_contributions()` [packages/treetime/src/optimize/eval.rs#L24](../../packages/treetime/src/optimize/eval.rs#L24) computes `exp_ev = (eigvals * branch_length).mapv(f64::exp)`

The Newton and Brent methods ([packages/treetime/src/optimize/method_newton.rs](../../packages/treetime/src/optimize/method_newton.rs), [packages/treetime/src/optimize/method_brent.rs](../../packages/treetime/src/optimize/method_brent.rs)) and the grid search all evaluate through this path.

## Impact

Negligible with one partition per run: each edge then has one contribution, and nothing is recomputed. The cost grows with the number of partitions that share a GTR model. The vector has one entry per alphabet state, so each recomputation is small compared with the per-site work.

## Proposed solution

Compute $\exp(\lambda t)$ once per distinct GTR model and branch length in each evaluation, and pass it to the per-contribution evaluation.
