# The branch-length objective omits the GTR rate `mu`

The per-branch objective of `optimize` and of the timetree branch-length likelihoods drops the overall rate of the GTR. It matches the tree likelihood only when `mu = 1`.

## Problem

`fn evaluate_site_contributions` [packages/treetime/src/optimize/eval.rs#L20](../../packages/treetime/src/optimize/eval.rs#L20) evaluates the first expression below. `fn GTR::expQt` [packages/treetime/src/gtr/gtr.rs#L129](../../packages/treetime/src/gtr/gtr.rs#L129) and the marginal passes use the second:

$$
\ell_{\text{objective}}(t) = \sum_k c_k\, e^{\lambda_k t}, \qquad \ell_{\text{tree}}(t) = \sum_k c_k\, e^{\mu \lambda_k t}
$$

where:

- $t$ -- branch length
- $c_k$ -- coefficient of eigenmode $k$, built from the two messages on the branch
- $\lambda_k$ -- eigenvalue $k$ of the normalized rate matrix (`gtr.eigvals`)
- $\mu$ -- overall rate of the GTR (`gtr.mu`)

Named models also differ from v0 in $\mu$. v1 multiplies the requested rate by the average rate of the unnormalized matrix [packages/treetime/src/gtr/gtr.rs#L62-L63](../../packages/treetime/src/gtr/gtr.rs#L62-L63), so JC69 gets $\mu = 0.75$. v0 `GTR.standard` sets $\mu$ back to the requested value, default 1, after it builds the model [packages/legacy/treetime/treetime/gtr.py#L460](../../packages/legacy/treetime/treetime/gtr.py#L460).

## Impact

- `optimize --model infer`: `fn normalize_partition_rates` [packages/treetime/src/optimize/run_loop.rs#L436](../../packages/treetime/src/optimize/run_loop.rs#L436) rescales to $\mu = 1$, so the objective is correct
- `optimize` with a named model: the two errors cancel, and the branch lengths equal v0's. The reported likelihood uses $\mu = 0.75$ for JC69
- `timetree`: the branch-length likelihoods [packages/treetime/src/timetree/inference/branch_length_likelihood.rs](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs) use the same contributions without rescaling, so their time scale is off by $\mu$ (about 0.9 to 1.1 for inferred models)

## Fix direction

Match v0. Both changes are needed together:

- Multiply the eigenvalues by $\mu$ in the objective: [packages/treetime/src/optimize/dense_eval.rs#L12](../../packages/treetime/src/optimize/dense_eval.rs#L12), [packages/treetime/src/optimize/sparse_eval.rs#L15](../../packages/treetime/src/optimize/sparse_eval.rs#L15)
- Give named models the requested $\mu$ (default 1) after construction, as v0 does. The constructors start at `fn jc69` [packages/treetime/src/gtr/get_gtr.rs#L122](../../packages/treetime/src/gtr/get_gtr.rs#L122)

The first change alone moves the total branch length of `test_gm_optimize` (JC69) from v0's 0.128 to 0.170 = 0.128 / 0.75.

## Validation

- `test_branch_objective_matches_tree_likelihood_with_skewed_pi` [packages/treetime/src/optimize/__tests__/test_objective_matches_tree_likelihood.rs#L25](../../packages/treetime/src/optimize/__tests__/test_objective_matches_tree_likelihood.rs#L25) passes without its `mu = 1` workaround and with a case where $\mu \neq 1$
- `test_gm_optimize` passes unchanged
