# Branch likelihood evaluation recomputes exp(eigvals\*t) per contribution

Every evaluation of the branch-length objective computes the vector $\exp(\lambda t)$ of the GTR eigenvalues $\lambda$ at branch length $t$ once for each `OptimizationContribution`. When several contributions on one edge share the same GTR eigenvalues, the same vector is computed several times per evaluation.

## Locations

- `fn sum_partition_metrics()` [packages/treetime/src/optimize/likelihood.rs#L42-L56](../../packages/treetime/src/optimize/likelihood.rs#L42-L56) loops over the contributions of an edge and evaluates each one separately
- `fn evaluate_dense_contribution()` and `fn evaluate_sparse_contribution()` pass the eigenvalues of the contribution's own GTR [packages/treetime/src/optimize/dense_eval.rs#L12](../../packages/treetime/src/optimize/dense_eval.rs#L12), [packages/treetime/src/optimize/sparse_eval.rs#L15](../../packages/treetime/src/optimize/sparse_eval.rs#L15)
- `fn evaluate_site_contributions()` [packages/treetime/src/optimize/eval.rs#L23](../../packages/treetime/src/optimize/eval.rs#L23) computes `exp_ev` from the rate eigenvalues of the contribution

The Newton and Brent methods ([packages/treetime/src/optimize/method_newton.rs](../../packages/treetime/src/optimize/method_newton.rs), [packages/treetime/src/optimize/method_brent.rs](../../packages/treetime/src/optimize/method_brent.rs)) and the grid search all evaluate through this path.

## Impact

Negligible with one partition per run: each edge then has one contribution, and nothing is recomputed. The cost grows with the number of partitions that share a GTR model. The vector has one entry per alphabet state, so each recomputation is small compared with the per-site work.

## Proposed solution

Compute $\exp(\lambda t)$ once per distinct GTR model and branch length in each evaluation, and pass it to the per-contribution evaluation.

## Proposed implementation plan

1. Profile first on a multi-partition run (`--partition` entries that share one model) to confirm that the per-contribution `exp_ev` in `fn evaluate_site_contributions()` (`eval.rs#L23`) shows in the profile; close the issue as negligible if it does not.
2. Split `evaluate_site_contributions()` into a step that builds `exp_ev` and `ev_exp_ev` from `rate_eigvals` and $t$, and a step that consumes them; in `fn sum_partition_metrics()` (`likelihood.rs#L42`) group the contributions by equal `rate_eigvals` and build the vectors once per group.
3. Test: an `rstest` that evaluates two contributions with equal eigenvalues through the grouped path and one by one, and asserts equal `OptimizationMetrics` with `assert_ulps_eq!` (oracle: the per-contribution evaluation, which computes the same terms independently per contribution).
4. Validate with `just lint-rs` and `just test-rs`; outputs are unchanged by construction.
