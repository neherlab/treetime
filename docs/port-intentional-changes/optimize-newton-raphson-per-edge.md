# Per-edge branch length optimization uses Newton-Raphson instead of Brent

v1 uses Newton-Raphson with analytical first and second derivatives for per-edge branch length optimization. v0 uses Brent's derivative-free method via `scipy.optimize.minimize_scalar`. Both operate on the same eigendecomposition-based likelihood function.

## What v0 does

v0 `optimal_t_compressed()` ([packages/legacy/treetime/treetime/gtr.py#L816-L920](../../packages/legacy/treetime/treetime/gtr.py#L816-L920)) optimizes each branch length using Brent's method (Brent, 1973) with three design choices:

1. **sqrt(t) reparameterization.** The objective function takes `s = sqrt(t)` as input and computes likelihood at `t = s^2`. This reshapes the likelihood surface near zero, where `dL/dt` can be steep. The bracket is `[-sqrt(MAX_BRANCH_LENGTH), sqrt(hamming_distance), sqrt(MAX_BRANCH_LENGTH)]`.

2. **Derivative-free.** Brent's method combines golden section search (linear convergence, ratio ~0.618) with parabolic interpolation (superlinear, order ~1.325). It needs only function evaluations, no derivatives.

3. **Hamming distance fallback.** If `scipy.optimize.minimize_scalar` reports failure, v0 falls back to the raw Hamming distance between parent and child sequences as the branch length estimate.

v0 also adds a regularization penalty `exp(t^4/10000)` when optimizing with marginal profiles (`profiles=True`) to prevent unbounded branch growth.

## What v1 does

v1 `run_optimize_mixed()` ([packages/treetime/src/commands/optimize/optimize_unified.rs](../../packages/treetime/src/commands/optimize/optimize_unified.rs)) offers three per-edge optimization methods selectable via `--opt-method`:

1. **`newton-sqrt` (default).** Newton-Raphson in $\sqrt{t}$ space. Reparameterizes as $s = \sqrt{t}$ and applies the chain rule to transform $t$-space derivatives. This reduces the indel Hessian singularity from $O(1/t^2)$ to $O(1/t)$, improving conditioning on short branches with indels. Falls back to grid search when the $s$-space Hessian is non-negative.

2. **`newton`.** Newton-Raphson in $t$ space with analytical derivatives from eigendecomposition-based coefficient caching. The Newton step is clamped to `[-1.0, current_bl]`. Max 10 inner iterations, convergence at `|delta_bl| <= max(0.001 * bl, 1e-8)`. Falls back to grid search when the Hessian is non-negative.

3. **`brent`.** Brent's method via `argmin::BrentOpt`. Derivative-free, bracket-based. The bracket lower bound is `min_branch_length` (or 1e-12 when no indels); the upper bound is `max(1.5 * bl + one_mutation, 0.5)`.

All methods share:

- **Grid search fallback** for non-concave regions (100 log-spaced points).
- **Zero-branch short-circuit** using the derivative-sign criterion at $t = 0$ for unimodal models (JC69, F81).

## Why v1 changes this

**Analytical derivatives are available at negligible cost.** v1's eigenvalue-space coefficient caching (`k_c`) computes the likelihood as `sum_c k_c exp(lambda_c * t)`. The first and second derivatives with respect to `t` follow by multiplying `k_c` by `lambda_c` and `lambda_c^2` respectively, reusing the same `exp(lambda_c * t)` values. Brent's method would discard these derivatives and perform 3-6x more function evaluations to reach the same precision.

**Newton-Raphson is the standard in phylogenetic ML software.** Both RAxML (Stamatakis, 2006; 2014) and IQ-TREE (Nguyen et al., 2015; Minh et al., 2020) use Newton-Raphson for per-branch optimization. IQ-TREE uses Brent as a fallback, matching v1's layered approach.

**Quadratic convergence.** Near the optimum, Newton-Raphson doubles the number of correct digits per iteration. Typical convergence: 3-5 iterations per branch. Brent achieves superlinear convergence of order ~1.325 (successive parabolic interpolation), requiring 10-30 function evaluations. For the ~2N branches in a tree with N leaves, this difference compounds.

**Unimodality for common models.** Dinh & Matsen (arXiv 1507.03647) prove the one-dimensional phylogenetic likelihood function has at most one stationary point under JC69, F81, and all binary models. For these models, Newton-Raphson converges to the global maximum from any starting point where the Hessian is negative. The grid search fallback handles the non-concave case that can arise under K80, HKY85, and GTR.

## Tradeoffs

**Recovered: sqrt(t) reparameterization.** The default `newton-sqrt` method restores v0's sqrt(t) conditioning. The chain rule transformation is $O(1)$ per evaluation.

**Still lost: regularization penalty.** v0 adds `exp(t^4/10000)` when optimizing with profiles. v1 does not apply this penalty. The grid search upper bound and Brent bracket provide a soft cap.

**Still lost: Hamming distance fallback.** v0 returns raw Hamming distance when Brent reports failure. v1 has no equivalent fallback: if optimization fails, the branch length remains at its current value.

**Gained: method selection.** Users can choose between three methods via `--opt-method`, matching the pattern in the `clock` command's `--branch-split-method`.

**Gained: cheaper per-evaluation cost (Newton methods).** Each Newton iteration gets likelihood AND derivatives from the same eigenvalue exponentials. Brent needs one full likelihood evaluation per bracket refinement step with no derivative reuse.

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
