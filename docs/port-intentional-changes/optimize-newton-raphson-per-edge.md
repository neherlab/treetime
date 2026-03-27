# Per-edge branch length optimization uses Newton-Raphson instead of Brent

v1 uses Newton-Raphson with analytical first and second derivatives for per-edge branch length optimization. v0 uses Brent's derivative-free method via `scipy.optimize.minimize_scalar`. Both operate on the same eigendecomposition-based likelihood function.

## What v0 does

v0 `optimal_t_compressed()` ([packages/legacy/treetime/treetime/gtr.py#L816-L920](../../packages/legacy/treetime/treetime/gtr.py#L816-L920)) optimizes each branch length using Brent's method (Brent, 1973) with three design choices:

1. **sqrt(t) reparameterization.** The objective function takes `s = sqrt(t)` as input and computes likelihood at `t = s^2`. This reshapes the likelihood surface near zero, where `dL/dt` can be steep. The bracket is `[-sqrt(MAX_BRANCH_LENGTH), sqrt(hamming_distance), sqrt(MAX_BRANCH_LENGTH)]`.

2. **Derivative-free.** Brent's method combines golden section search (linear convergence, ratio ~0.618) with parabolic interpolation (superlinear, order ~1.325). It needs only function evaluations, no derivatives.

3. **Hamming distance fallback.** If `scipy.optimize.minimize_scalar` reports failure, v0 falls back to the raw Hamming distance between parent and child sequences as the branch length estimate.

v0 also adds a regularization penalty `exp(t^4/10000)` when optimizing with marginal profiles (`profiles=True`) to prevent unbounded branch growth.

## What v1 does

v1 `run_optimize_mixed()` ([packages/treetime/src/commands/optimize/optimize_unified.rs#L221-L291](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L221-L291)) uses Newton-Raphson with two fallback layers:

1. **Newton-Raphson with analytical derivatives.** The eigendecomposition `Q = V diag(lambda) V^{-1}` allows computing likelihood, first derivative, and second derivative from the same cached coefficients `k_c = (msg.dot(V)) * (msg.dot(V_inv.T))` (dense: [packages/treetime/src/commands/optimize/optimize_dense.rs#L67-L75](../../packages/treetime/src/commands/optimize/optimize_dense.rs#L67-L75), sparse: [packages/treetime/src/commands/optimize/optimize_sparse.rs#L46-L119](../../packages/treetime/src/commands/optimize/optimize_sparse.rs#L46-L119)). The derivatives come at negligible cost beyond the likelihood evaluation:

   ```
   L   = sum_j log(sum_c k_c exp(lambda_c * t))
   L'  = sum_j (sum_c k_c lambda_c exp(lambda_c * t)) / (sum_c k_c exp(lambda_c * t))
   L'' = sum_j (sum_c k_c lambda_c^2 exp(lambda_c * t)) / (...) - (L')^2
   ```

   The Newton step is clamped to `[-1.0, current_bl]` to enforce non-negativity ([optimize_unified.rs#L254](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L254)). Max 10 inner iterations, convergence at `|delta_bl| <= 0.001 * bl`.

2. **Grid search fallback.** When the second derivative is non-negative (non-concave region where Newton would step the wrong way), v1 evaluates 100 equally-spaced points on `[0.1 * one_mutation, 1.5 * bl + one_mutation]` and selects the maximum ([optimize_unified.rs#L273-L286](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L273-L286)).

3. **Zero-branch short-circuit.** Before any optimization, v1 checks if zero branch length is optimal: combined likelihood at zero > 0.01 AND first derivative < 0 at zero (`is_zero_branch_optimal()` at [optimize_unified.rs#L192-L214](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L192-L214)).

## Why v1 changes this

**Analytical derivatives are available at negligible cost.** v1's eigenvalue-space coefficient caching (`k_c`) computes the likelihood as `sum_c k_c exp(lambda_c * t)`. The first and second derivatives with respect to `t` follow by multiplying `k_c` by `lambda_c` and `lambda_c^2` respectively, reusing the same `exp(lambda_c * t)` values. Brent's method would discard these derivatives and perform 3-6x more function evaluations to reach the same precision.

**Newton-Raphson is the standard in phylogenetic ML software.** Both RAxML (Stamatakis, 2006; 2014) and IQ-TREE (Nguyen et al., 2015; Minh et al., 2020) use Newton-Raphson for per-branch optimization. IQ-TREE uses Brent as a fallback, matching v1's layered approach.

**Quadratic convergence.** Near the optimum, Newton-Raphson doubles the number of correct digits per iteration. Typical convergence: 3-5 iterations per branch. Brent achieves superlinear convergence of order ~1.325 (successive parabolic interpolation), requiring 10-30 function evaluations. For the ~2N branches in a tree with N leaves, this difference compounds.

**Unimodality for common models.** Dinh & Matsen (arXiv 1507.03647) prove the one-dimensional phylogenetic likelihood function has at most one stationary point under JC69, F81, and all binary models. For these models, Newton-Raphson converges to the global maximum from any starting point where the Hessian is negative. The grid search fallback handles the non-concave case that can arise under K80, HKY85, and GTR.

## Tradeoffs

**Lost: sqrt(t) reparameterization.** v0's optimization in sqrt(t) space improves conditioning near zero where dL/dt can be large. v1 compensates with the zero-branch short-circuit and Newton step clamping, but does not reshape the likelihood surface itself.

**Lost: regularization penalty.** v0 adds `exp(t^4/10000)` when optimizing with profiles. v1 does not apply this penalty. For well-behaved datasets this is immaterial; for cases where the likelihood surface is flat at large branch lengths, v0's penalty prevents runaway. The grid search upper bound `1.5 * bl + one_mutation` provides a soft cap in v1.

**Lost: Hamming distance fallback.** v0 returns raw Hamming distance when Brent reports failure. v1 has no equivalent fallback: if both Newton and grid search fail to improve the likelihood, the branch length remains at its current value.

**Gained: cheaper per-evaluation cost.** Each Newton iteration gets likelihood AND derivatives from the same eigenvalue exponentials. Brent needs one full likelihood evaluation per bracket refinement step with no derivative reuse.

## Practical impact

For end-users, the per-edge optimization method is invisible. Branch lengths converge to the same ML estimates (within numerical tolerance) regardless of whether Newton or Brent finds them, given the same GTR model and marginal profiles.

The convergence issues previously reported for the optimize command (oscillation between iterations) stemmed from the undamped outer loop, not from the per-branch optimizer. v1 now applies the same exponential damping as v0 (`--damping`, default 0.75). Brent per-branch in v0 is orthogonal to this.

## References

- Brent RP (1973). _Algorithms for Minimization without Derivatives_. Prentice-Hall.
- Dinh V, Matsen FA IV (2015). The shape of the one-dimensional phylogenetic likelihood function. arXiv 1507.03647.
- Stamatakis A (2006). RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models. Bioinformatics 22(21):2688-2690.
- Stamatakis A (2014). RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics 30(9):1312-1313.
- Nguyen LT, Schmidt HA, von Haeseler A, Minh BQ (2015). IQ-TREE: A fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Mol Biol Evol 32(1):268-274.
- Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Mol Biol Evol 37(5):1530-1534.
