# Per-edge branch length optimization: Brent and Newton-Raphson in sqrt(t) space

v1 offers two per-edge branch length optimization methods via `--opt-method`: Brent's method and Newton-Raphson, both in $\sqrt{t}$ space. The default `brent` uses the v0 method, Brent in $\sqrt{t}$ space, but differs from v0 in the search interval, the start, the objective penalty, the stopping tolerance, and the failure handling, so per-edge results differ slightly ([kb/issues/M-optimize-brent-setup-differs-from-v0.md](../issues/M-optimize-brent-setup-differs-from-v0.md)). v0 offers only Brent in $\sqrt{t}$ space.

## What v0 does

v0 `optimal_t_compressed()` ([packages/legacy/treetime/treetime/gtr.py#L816-L920](../../packages/legacy/treetime/treetime/gtr.py#L816-L920)) optimizes each branch length using Brent's method (Brent, 1973) with three design choices:

1. sqrt(t) reparameterization. The objective function takes `s = sqrt(t)` as input and computes likelihood at `t = s^2`. This reshapes the likelihood surface near zero, where `dL/dt` can be steep. The bracket is `[-sqrt(MAX_BRANCH_LENGTH), sqrt(hamming_distance), sqrt(MAX_BRANCH_LENGTH)]`.

2. Derivative-free. Brent's method combines golden section search (linear convergence, ratio ~0.618) with parabolic interpolation (superlinear, order ~1.325). It needs only function evaluations, no derivatives.

3. Hamming distance fallback. If `scipy.optimize.minimize_scalar` reports failure, v0 falls back to the raw Hamming distance between parent and child sequences as the branch length estimate.

v0 also adds a regularization penalty `exp(s^4/10000)`, that is $e^{t^2/10^4}$, when optimizing with marginal profiles (`profiles=True`) to prevent unbounded branch growth. In the code the argument of `_neg_prob` is named `t` but holds $s$.

## What v1 does

v1 `run_optimize_mixed()` ([packages/treetime/src/optimize/dispatch.rs](../../packages/treetime/src/optimize/dispatch.rs)) offers two per-edge optimization methods selectable via `--opt-method`:

### Brent's method (derivative-free)

1. `brent` (default). Brent's method in $\sqrt{t}$ space via `argmin::BrentOpt`. The cost function evaluates $-\ell(s^2)$ and the interval endpoints transform as $\sqrt{\text{lower}}$, $\sqrt{\text{upper}}$, as in v0. It differs from v0 on the success path too: the interval is $[\max(t_\text{min}, 10^{-12}), \max(1.5\,t_0 + \text{one mutation}, 0.5)]$ instead of $t \le 4$, the solver starts at a golden-section point instead of the Hamming distance and can return a point slightly worse than the start, the objective has no penalty, and the tolerance is absolute ($10^{-5}$ in $s$) instead of v0's relative per-iteration schedule ([kb/issues/M-optimize-brent-setup-differs-from-v0.md](../issues/M-optimize-brent-setup-differs-from-v0.md)). Tolerance $\epsilon_s$ in $s$-space maps to $t$-space precision $\approx 2s^* \epsilon_s$, tighter near zero. Differs from v0 on solver failure: v1 surfaces an error to the caller via `Result<f64, Report>`, while v0 silently substitutes the raw Hamming distance when SciPy reports non-convergence. v1 does not read argmin's termination status, so a stop at the iteration limit is used silently.

### Newton-Raphson (analytical derivatives)

The substitution Hessian (posterior variance of eigenvalues) is computed in the centered Welford form $\sum_c w_c (\lambda_c - \bar\lambda)^2$, which avoids the catastrophic cancellation of the two-moment form $E[\lambda^2] - E[\lambda]^2$ when the posterior is tightly peaked.

2. `newton`. Newton-Raphson in $\sqrt{t}$ space with analytical derivatives from eigendecomposition-based coefficient caching. Reparameterizes as $s = \sqrt{t}$ and applies the chain rule: $d\ell/ds = 2s \cdot d\ell/dt$, $d^2\ell/ds^2 = 4s^2 \cdot d^2\ell/dt^2 + 2 \cdot d\ell/dt$. Reduces the Poisson indel Hessian singularity from $O(1/t^2)$ in $t$ space to $O(1/t)$, so the indel term dominates the substitution Hessian on short branches less than in $t$ space.

Both methods share:

- Grid search fallback for non-concave regions (Newton-Raphson only, 100 log-spaced points).
- Zero-branch short-circuit using the derivative-sign criterion at $t = 0$ for unimodal models (JC69, F81).

## Why v1 offers method selection

**One parameterization.** The $\sqrt{t}$ parameterization is the v0 parameterization and improves numerical conditioning on short branches and indel-bearing edges. Both methods use it, so a change of `--opt-method` changes only the solver.

**Newton-Raphson for analytical derivatives.** v1's eigenvalue-space coefficient caching (`k_c`) computes the likelihood as `sum_c k_c exp(lambda_c * t)`. The first and second derivatives follow by multiplying `k_c` by `lambda_c` and `lambda_c^2`, reusing the same exponentials. Newton achieves quadratic convergence near the optimum (3-5 iterations typical vs 10-30 function evaluations for Brent).

**Brent for robustness.** Brent's method is derivative-free and guaranteed to converge within the bracket. It serves as both a reference implementation and a fallback for cases where Newton behaves poorly.

**Default is the v0 method.** The `brent` default keeps v0's solver and parameterization without configuration; the setup differences from v0 are recorded in [kb/issues/M-optimize-brent-setup-differs-from-v0.md](../issues/M-optimize-brent-setup-differs-from-v0.md).

## Tradeoffs

**sqrt(t) reparameterization.** Used by `brent` (default; the interval, start, objective, tolerance, and failure handling differ from v0, see [kb/issues/M-optimize-brent-setup-differs-from-v0.md](../issues/M-optimize-brent-setup-differs-from-v0.md)) and `newton`.

**Still lost: regularization penalty.** v0 adds `exp(s^4/10000)` ($e^{t^2/10^4}$) when optimizing with profiles. v1 does not apply this penalty. The grid search upper bound and Brent bracket provide a soft cap.

**Still lost: Hamming distance fallback.** v0 returns raw Hamming distance when Brent reports failure. v1 has no equivalent: solver errors propagate through `Result<f64, Report>` to `run_optimize_mixed`, which surfaces them to the caller; a stop at the iteration limit is not detected. v1 prefers a loud error over v0's silent substitution because the silent fallback biased downstream branch lengths and made likelihood comparisons unreliable.

**Gained: method selection.** Users can choose between Brent's method and Newton-Raphson via `--opt-method`, matching the pattern in the `clock` command's `--branch-split-method`.

## Practical impact

For end-users, the per-edge optimization method is invisible. Branch lengths converge to the same ML estimates (within numerical tolerance) regardless of method, given the same GTR model and marginal profiles. Both methods produce log-likelihood values agreeing within 1e-3 on indel-bearing edges.

The outer loop applies no damping by default (`--damping`, default 0; v0 damps with 0.75). Optimizing edges in series as one sum removes the 2-cycle that an undamped loop with independent edge updates would produce ([optimize-damping-disabled-by-default.md](optimize-damping-disabled-by-default.md), [optimize-series-edges-joint-sum.md](optimize-series-edges-joint-sum.md)).

## References

- Brent, Richard P. 1973. _Algorithms for Minimization Without Derivatives._ Prentice-Hall. ISBN 978-0-13-022335-7.
- Dinh, Vu, and Frederick A. Matsen IV. 2015. "The Shape of the One-Dimensional Phylogenetic Likelihood Function." arXiv:1507.03647. https://arxiv.org/abs/1507.03647
- Felsenstein, Joseph. 2003. _Inferring Phylogenies._ Sinauer Associates. ISBN 978-0-87893-177-4. Chapter 16.
- Stamatakis, Alexandros. 2006. "RAxML-VI-HPC: Maximum Likelihood-Based Phylogenetic Analyses with Thousands of Taxa and Mixed Models." _Bioinformatics_ 22(21):2688-2690. https://doi.org/10.1093/bioinformatics/btl446
- Stamatakis, Alexandros. 2014. "RAxML Version 8: A Tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies." _Bioinformatics_ 30(9):1312-1313. https://doi.org/10.1093/bioinformatics/btu033
- Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh. 2015. "IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies." _Molecular Biology and Evolution_ 32(1):268-274. https://doi.org/10.1093/molbev/msu300
- Minh, Bui Quang, Heiko A. Schmidt, Olga Chernomor, Dominik Schrempf, Michael D. Woodhams, Arndt von Haeseler, and Robert Lanfear. 2020. "IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era." _Molecular Biology and Evolution_ 37(5):1530-1534. https://doi.org/10.1093/molbev/msaa015
