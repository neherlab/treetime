# The branch-length objective ignores the substitution rate of the GTR

`fn evaluate_site_contributions` in [packages/treetime/src/optimize/eval.rs](../../packages/treetime/src/optimize/eval.rs) evaluates each site as $\sum_k c_k e^{\lambda_k t}$, while the marginal passes use `GTR::expQt`, $e^{\mu \lambda_k t}$. The per-branch objective therefore matches the tree likelihood only when `gtr.mu == 1`.

## Where it matters

- `optimize` rescales the model to `mu = 1` (`fn normalize_partition_rates`) only for `--model infer`. With a named model it keeps the model's `mu`: JC69 in v1 has `mu = 0.75` (average substitution rate 0.75), JC69 in v0 has `mu = 1.0`. The omission makes `optimize --model jc69` return v0's branch lengths, while its reported likelihood uses `mu = 0.75`.
- `timetree` builds branch-length likelihoods from the same contributions ([packages/treetime/src/timetree/inference/branch_length_likelihood.rs](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs)) and does not normalize the rate, so they are on a time scale off by `mu` (inferred models have `mu` around 0.9 to 1.1).

## Evidence

Multiplying the eigenvalues by `mu` in `dense_eval.rs` and `sparse_eval.rs` makes the objective match the tree likelihood exactly (`test_branch_objective_matches_tree_likelihood_with_skewed_pi`, which currently sets `mu = 1`), but moves `test_gm_optimize` (JC69, no normalization) from v0's total branch length 0.128 to 0.170 = 0.128 / 0.75.

## Fix direction

- Evaluate each site as $\sum_k c_k e^{\mu \lambda_k t}$, as `GTR::expQt` does, so that the objective matches the tree likelihood for any `mu`
- Give named models the requested `mu` (default 1) after construction, as v0 does: `GTR.standard` builds the model and then sets `model.mu` ([packages/legacy/treetime/treetime/gtr.py#L460](../../packages/legacy/treetime/treetime/gtr.py#L460)). The v1 `GTR` builder keeps `mu` multiplied by the average rate of `W` ([packages/treetime/src/gtr/gtr.rs#L63](../../packages/treetime/src/gtr/gtr.rs#L63)), which gives JC69 `mu = 0.75`
- With both changes, `test_gm_optimize` stays at v0's total branch length, and `test_branch_objective_matches_tree_likelihood_with_skewed_pi` no longer needs to set `mu = 1`
