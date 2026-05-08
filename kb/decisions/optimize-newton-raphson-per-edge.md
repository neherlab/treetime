# Per-edge branch length optimization: 6-method selection

v1 offers six per-edge branch length optimization methods via `--opt-method`, combining two algorithms (Newton-Raphson, Brent) with three parameterizations ($t$, $\sqrt{t}$, $\ln(t)$). The default `brent-sqrt` matches v0 exactly. v0 offers only Brent in $\sqrt{t}$ space.

## What v0 does

v0 `optimal_t_compressed()` ([packages/legacy/treetime/treetime/gtr.py#L816-L920](../../packages/legacy/treetime/treetime/gtr.py#L816-L920)) optimizes each branch length using Brent's method (Brent, 1973) with three design choices:

1. sqrt(t) reparameterization. The objective function takes `s = sqrt(t)` as input and computes likelihood at `t = s^2`. This reshapes the likelihood surface near zero, where `dL/dt` can be steep. The bracket is `[-sqrt(MAX_BRANCH_LENGTH), sqrt(hamming_distance), sqrt(MAX_BRANCH_LENGTH)]`.

2. Derivative-free. Brent's method combines golden section search (linear convergence, ratio ~0.618) with parabolic interpolation (superlinear, order ~1.325). It needs only function evaluations, no derivatives.

3. Hamming distance fallback. If `scipy.optimize.minimize_scalar` reports failure, v0 falls back to the raw Hamming distance between parent and child sequences as the branch length estimate.

v0 also adds a regularization penalty `exp(t^4/10000)` when optimizing with marginal profiles (`profiles=True`) to prevent unbounded branch growth.

## What v1 does

v1 `run_optimize_mixed()` ([packages/treetime/src/commands/optimize/optimize_unified.rs](../../packages/treetime/src/commands/optimize/optimize_unified.rs)) offers six per-edge optimization methods selectable via `--opt-method`:

### Brent methods (derivative-free)

1. `brent-sqrt` (default). Brent's method in $\sqrt{t}$ space via `argmin::BrentOpt`. Matches v0 on the success path: the cost function evaluates $-\ell(s^2)$ and bracket endpoints transform as $\sqrt{\text{lower}}$, $\sqrt{\text{upper}}$. Tolerance $\epsilon_s$ in $s$-space maps to $t$-space precision $\approx 2s^* \epsilon_s$, tighter near zero. Differs from v0 on solver failure: v1 surfaces the error to the caller via `Result<f64, Report>`, while v0 silently substitutes the raw Hamming distance.

2. `brent`. Brent's method in $t$ space. Derivative-free, bracket-based. Included for completeness; `brent-sqrt` dominates for convergence speed near $t = 0$.

3. `brent-log`. Brent's method in $\ln(t)$ space. The cost function evaluates $-\ell(e^u)$ and bracket endpoints transform as $\ln(\text{lower})$, $\ln(\text{upper})$. Tolerance $\epsilon_u$ in $u$-space is a natural relative tolerance ($dt/t \approx du$). Smoothest objective surface of all Brent variants.

### Newton methods (analytical derivatives)

The substitution Hessian (posterior variance of eigenvalues) is computed in the centered Welford form $\sum_c w_c (\lambda_c - \bar\lambda)^2$, which avoids the catastrophic cancellation of the two-moment form $E[\lambda^2] - E[\lambda]^2$ when the posterior is tightly peaked.

4. `newton`. Newton-Raphson in $t$ space with analytical derivatives from eigendecomposition-based coefficient caching. Baseline matching RAxML-NG/IQ-TREE. The Poisson indel Hessian ($-k/t^2$) can dominate on short branches with indels.

5. `newton-sqrt`. Newton-Raphson in $\sqrt{t}$ space. Reparameterizes as $s = \sqrt{t}$ and applies the chain rule: $d\ell/ds = 2s \cdot d\ell/dt$, $d^2\ell/ds^2 = 4s^2 \cdot d^2\ell/dt^2 + 2 \cdot d\ell/dt$. Reduces the indel Hessian singularity from $O(1/t^2)$ to $O(1/t)$.

6. `newton-log`. Newton-Raphson in $\ln(t)$ space. Reparameterizes as $u = \ln(t)$ and applies the chain rule: $d\ell/du = t \cdot d\ell/dt$, $d^2\ell/du^2 = t^2 \cdot d^2\ell/dt^2 + t \cdot d\ell/dt$. Eliminates the indel Hessian singularity entirely ($\ell''_{\text{indel}} = -\mu t$, bounded). Best conditioning of all Newton variants.

All methods share:

- Grid search fallback for non-concave regions (Newton methods only, 100 log-spaced points).
- Zero-branch short-circuit using the derivative-sign criterion at $t = 0$ for unimodal models (JC69, F81).

## Why v1 offers method selection

**Different tradeoffs for different use cases.** The $\sqrt{t}$ and $\ln(t)$ parameterizations improve numerical conditioning on short branches and indel-bearing edges, but add transformation overhead. Pure $t$-space methods are simpler and sufficient for well-conditioned cases.

**Newton-Raphson for analytical derivatives.** v1's eigenvalue-space coefficient caching (`k_c`) computes the likelihood as `sum_c k_c exp(lambda_c * t)`. The first and second derivatives follow by multiplying `k_c` by `lambda_c` and `lambda_c^2`, reusing the same exponentials. Newton achieves quadratic convergence near the optimum (3-5 iterations typical vs 10-30 function evaluations for Brent).

**Brent for robustness.** Brent's method is derivative-free and guaranteed to converge within the bracket. It serves as both a reference implementation and a fallback for cases where Newton behaves poorly.

**Default matches v0.** The `brent-sqrt` default ensures golden master tests pass and users get v0-equivalent behavior without configuration.

## Tradeoffs

**sqrt(t) reparameterization.** Available as `brent-sqrt` (default, matches v0 on the success path; on solver failure v1 returns `Err(Report)` rather than v0's silent Hamming-distance fallback) and `newton-sqrt`.

**New: ln(t) reparameterization.** Available as `brent-log` and `newton-log`. Eliminates the indel Hessian singularity entirely. Not present in v0 or other phylogenetic tools for branch length optimization. Precedented by coalescent Tc optimization in this codebase ([packages/treetime/src/commands/timetree/coalescent/optimize_tc.rs#L62](../../packages/treetime/src/commands/timetree/coalescent/optimize_tc.rs#L62)).

**Still lost: regularization penalty.** v0 adds `exp(t^4/10000)` when optimizing with profiles. v1 does not apply this penalty. The grid search upper bound and Brent bracket provide a soft cap.

**Still lost: Hamming distance fallback.** v0 returns raw Hamming distance when Brent reports failure. v1 has no equivalent: solver failures propagate through `Result<f64, Report>` to `run_optimize_mixed`, which surfaces them to the caller. v1 prefers a loud error over v0's silent substitution because the silent fallback biased downstream branch lengths and made likelihood comparisons unreliable.

**Gained: method selection.** Users can choose between six methods via `--opt-method`, matching the pattern in the `clock` command's `--branch-split-method`.

## Practical impact

For end-users, the per-edge optimization method is invisible. Branch lengths converge to the same ML estimates (within numerical tolerance) regardless of method, given the same GTR model and marginal profiles. All six methods produce log-likelihood values agreeing within 1e-3 on indel-bearing edges.

Convergence behavior tracks v0: the outer loop applies exponential damping (`--damping`, default 0.75) so alternating reconstruction and per-branch optimization stay bounded against the 2-cycle that an undamped loop would produce.

## References

- Brent, Richard P. 1973. _Algorithms for Minimization Without Derivatives._ Prentice-Hall. ISBN 978-0-13-022335-7.
- Dinh, Vu, and Frederick A. Matsen IV. 2015. "The Shape of the One-Dimensional Phylogenetic Likelihood Function." arXiv:1507.03647. https://arxiv.org/abs/1507.03647
- Felsenstein, Joseph. 2003. _Inferring Phylogenies._ Sinauer Associates. ISBN 978-0-87893-177-4. Chapter 16.
- Stamatakis, Alexandros. 2006. "RAxML-VI-HPC: Maximum Likelihood-Based Phylogenetic Analyses with Thousands of Taxa and Mixed Models." _Bioinformatics_ 22(21):2688-2690. https://doi.org/10.1093/bioinformatics/btl446
- Stamatakis, Alexandros. 2014. "RAxML Version 8: A Tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies." _Bioinformatics_ 30(9):1312-1313. https://doi.org/10.1093/bioinformatics/btu033
- Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh. 2015. "IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies." _Molecular Biology and Evolution_ 32(1):268-274. https://doi.org/10.1093/molbev/msu300
- Minh, Bui Quang, Heiko A. Schmidt, Olga Chernomor, Dominik Schrempf, Michael D. Woodhams, Arndt von Haeseler, and Robert Lanfear. 2020. "IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era." _Molecular Biology and Evolution_ 37(5):1530-1534. https://doi.org/10.1093/molbev/msaa015
