# Chapter 2: Branch length optimization

[Back to index](_index.md) | Previous: [Chapter 1: Introduction](1-introduction.md) | Next: [Chapter 3: Coalescent skyline](3-coalescent-skyline.md)

The per-edge branch length optimization is the innermost and most frequently called optimization in the pipeline. Given fixed ancestral state distributions at both endpoints of an edge, find the branch length `t` that maximizes the likelihood. This is a 1D scalar optimization problem. See [Iterative tree refinement, Chapter 5](../iterative-tree-refinement/5-branch-length-optimization.md) for the mathematical formulation.

## Tool comparison

| Property                 | RAxML-NG                      | IQ-TREE 2                  | PhyML                | TreeTime v0            | TreeTime v1           |
| :----------------------- | :---------------------------- | :------------------------- | :------------------- | :--------------------- | :-------------------- |
| **Method**               | Newton-Raphson                | Newton-Raphson + bisection | Cubic Hermite spline | Brent (scipy)          | Newton-Raphson + grid |
| **Derivatives**          | 1st + 2nd analytical          | 1st + 2nd analytical       | 1st only analytical  | None (derivative-free) | 1st + 2nd analytical  |
| **Non-concave fallback** | `dx = -f/\|df\|` (safe NR)    | Bisection within bracket   | None (asserts)       | N/A (bracket-based)    | Grid search (100 pts) |
| **Max inner iterations** | 30                            | 100 (10 during NNI)        | 1000+20              | scipy default          | 10                    |
| **Branch bounds**        | [1e-4, 100]                   | [1e-6, 10]                 | [1e-8, 100]          | [0, 4] (sqrt space)    | [0, inf)              |
| **Convergence**          | \|f\| < 1e-4 or \|dx\| < 1e-4 | \|dx\| < min_bl            | \|dLH\| < 1e-3       | tol=1e-10              | \|dt\| < 0.001\*t     |
| **Outer smoothings**     | Up to 32                      | Up to 100                  | Until dLH < tol      | EM loop (max_iter)     | EM loop (max_iter)    |
| **Damping**              | None (SAFE mode reverts)      | Rollback if LH drops       | Tracks best LH       | Inline per-edge        | Exponential decay     |
| **Precomputed per edge** | sumtable (eigenbasis)         | theta buffer               | dot_prod array       | N/A                    | coefficients matrix   |
| **Reparameterization**   | None                          | None                       | None                 | sqrt(t)                | None                  |

## Newton-Raphson

The standard method in the field. The update step:

```
t_new = t - f'(t) / f''(t)
```

where `f'` and `f''` are the first and second derivatives of the log-likelihood with respect to branch length. Both are computed analytically from the eigendecomposition of the rate matrix (see [Iterative tree refinement, Chapter 2](../iterative-tree-refinement/2-substitution-models.md)). The key insight: profile-eigenvector dot products are branch-length-independent and computed once per edge, so each Newton iteration requires only O(s) multiplications per site (where s is the number of states).

Newton-Raphson converges quadratically near the optimum: each iteration doubles the number of correct digits. Typical convergence in 3-8 iterations (5-16 likelihood-equivalents including derivative cost). This is why RAxML-NG, IQ-TREE, and PAML all chose NR as the default (<a id="cite-1"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[1](#ref-1)]; <a id="cite-2"></a>[Stamatakis 2014](https://doi.org/10.1093/bioinformatics/btu033) [[2](#ref-2)]; <a id="cite-3"></a>[Nguyen et al. 2015](https://doi.org/10.1093/molbev/msu300) [[3](#ref-3)]).

**Non-concave fallback**: When `f''(t) >= 0`, the surface is not concave and Newton would step toward a minimum. Each tool handles this differently:

- **RAxML-NG** (coraxlib): uses `dx = -f / |df|`, preserving the NR direction but preventing sign inversion. Also has SAFE mode that recomputes likelihood and reverts if it decreased.
- **IQ-TREE**: falls back to bisection within the sign-change bracket `[xl, xh]`. Adapted from Numerical Recipes `rtsafe`.
- **TreeTime v1**: falls back to grid search (100 logarithmically spaced points via `geomspace`). This is the least principled fallback - BrentOpt would be more efficient (see [audit proposal P5](7-audit.md)).

## Brent's method

A derivative-free bracket-based optimizer combining golden-section search with successive parabolic interpolation (<a id="cite-4"></a>[Brent 1973](https://maths-people.anu.edu.au/~brent/pd/rpb011i.pdf) [[4](#ref-4)]). Convergence order ~1.325 (plastic number), falling back to golden-section ratio 0.618 on difficult surfaces.

TreeTime v0 is the only major phylogenetic tool using Brent for branch lengths, via `scipy.optimize.minimize_scalar(method='bounded')`. v0 adds a sqrt(t) reparameterization that spreads out the near-zero region and compresses the large-t region.

Brent requires 10-30 function evaluations (vs 3-8 NR iterations at ~2x cost each), so NR is typically faster when derivatives are available. Brent is preferred when derivatives are expensive or unavailable.

## Cubic Hermite spline (PhyML)

PhyML uses a cubic Hermite spline interpolation on the first derivative, replacing the legacy Brent implementation. Given two bracketing points `(u, fu, dfu)` and `(v, fv, dfv)` with function values and first derivatives, a cubic Hermite polynomial is constructed and its stationary points found via the quadratic formula.

This achieves cubic convergence using only first derivatives at two points - no second derivatives needed. The trade-off: maintaining two evaluation points per iteration instead of one. PhyML computes only first derivatives (via the eigendecomposition), not second derivatives. Source: `src/optimiz.c:2243` `Br_Len_Spline()` in [PhyML](https://github.com/stephaneguindon/phyml).

The method was described as "Newton-Raphson with cubic interpolation" in <a id="cite-5"></a>[Guindon and Gascuel 2003](https://doi.org/10.1080/10635150390235520) [[5](#ref-5)], but the actual implementation differs: it is not Newton-Raphson (which requires second derivatives at a single point), and has no graceful fallback - negative discriminants cause `assert(FALSE)`.

## Unimodality

The choice of optimizer matters because the 1D likelihood surface can be multimodal under some substitution models.

<a id="cite-6"></a>[Dinh and Matsen 2017](https://doi.org/10.1214/16-AAP1240) [[6](#ref-6)] proved:

- **JC69, F81, binary symmetric**: the 1D likelihood has at most one stationary point. Any optimizer converges to the global MLE.
- **K2P and above (HKY, GTR)**: explicit counterexamples with two local maxima exist. The space of rescaled 1D likelihoods under K2P is dense in all non-negative continuous functions on [0, infinity). Any shape is possible.

For TreeTime's common use case (JC69 on viral data), NR is provably safe. For GTR models, a grid-search or multi-start fallback is warranted. <a id="cite-7"></a>[Claywell et al. 2017](https://doi.org/10.1093/molbev/msx253) [[7](#ref-7)] developed a four-parameter surrogate function for 1D phylogenetic likelihoods that provides a theoretically grounded fast approximation.

## Code locations

**RAxML-NG** (via [coraxlib](https://codeberg.org/Exelixis-Lab/coraxlib) / [pll-modules](https://github.com/ddarriba/pll-modules)):

- Newton solver: `newton.c` `corax_opt_minimize_newton_multi()`
- Derivative kernel: `core_derivatives.c` (SIMD-vectorized sumtable + diagptable)
- Branch optimization driver: `opt_branches.c` `recomp_iterative_multi()`
- RAxML-NG wrapper: [src/TreeInfo.cpp#L339-L400](https://github.com/amkozlov/raxml-ng/blob/master/src/TreeInfo.cpp#L339-L400)

**IQ-TREE 2** ([iqtree2](https://github.com/iqtree/iqtree2)):

- Newton+bisection solver: [utils/optimization.cpp#L421-L498](https://github.com/iqtree/iqtree2/blob/master/utils/optimization.cpp#L421-L498) `minimizeNewton()`
- Derivative computation: [tree/phylokernelnew.h#L2235-L2550](https://github.com/iqtree/iqtree2/blob/master/tree/phylokernelnew.h#L2235-L2550) `computeLikelihoodDervSIMD()`
- Per-branch entry: [tree/phylotree.cpp#L2593-L2643](https://github.com/iqtree/iqtree2/blob/master/tree/phylotree.cpp#L2593-L2643) `optimizeOneBranch()`

**PhyML** ([phyml](https://github.com/stephaneguindon/phyml)):

- Cubic spline solver: [src/optimiz.c#L2243](https://github.com/stephaneguindon/phyml/blob/master/src/optimiz.c#L2243) `Br_Len_Spline()`
- Derivative: [src/lk.c#L655](https://github.com/stephaneguindon/phyml/blob/master/src/lk.c#L655) `dLk()`

**TreeTime v0** ([treetime](https://github.com/neherlab/treetime)):

- Branch optimization: [treetime/gtr.py#L816-L920](https://github.com/neherlab/treetime/blob/master/treetime/gtr.py#L816-L920) `optimal_t_compressed()`

**TreeTime v1**:

- Newton-Raphson + grid: [packages/treetime/src/commands/optimize/optimize_unified.rs#L221-L291](../../../packages/treetime/src/commands/optimize/optimize_unified.rs#L221-L291) `run_optimize_mixed()`

## References

1. <a id="ref-1"></a> Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _J. Mol. Evol._ 17:368-376. https://doi.org/10.1007/BF01734359 [↩](#cite-1)
2. <a id="ref-2"></a> Stamatakis, Alexandros. 2014. "RAxML Version 8." _Bioinformatics_ 30(9):1312-1313. https://doi.org/10.1093/bioinformatics/btu033 [↩](#cite-2)
3. <a id="ref-3"></a> Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh. 2015. "IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies." _Mol. Biol. Evol._ 32(1):268-274. https://doi.org/10.1093/molbev/msu300 [↩](#cite-3)
4. <a id="ref-4"></a> Brent, Richard P. 1973. _Algorithms for Minimization without Derivatives._ Prentice-Hall. https://maths-people.anu.edu.au/~brent/pd/rpb011i.pdf [↩](#cite-4)
5. <a id="ref-5"></a> Guindon, Stephane, and Olivier Gascuel. 2003. "A Simple, Fast, and Accurate Algorithm to Estimate Large Phylogenies by Maximum Likelihood." _Syst. Biol._ 52(5):696-704. https://doi.org/10.1080/10635150390235520 [↩](#cite-5)
6. <a id="ref-6"></a> Dinh, Vu, and Frederick A. Matsen IV. 2017. "The Shape of the One-Dimensional Phylogenetic Likelihood Function." _Ann. Appl. Prob._ 27(3):1646-1677. https://doi.org/10.1214/16-AAP1240 [↩](#cite-6)
7. <a id="ref-7"></a> Claywell, Brian C., Vu Dinh, Mathieu Fourment, Connor O. McCoy, and Frederick A. Matsen IV. 2017. "A Surrogate Function for One-Dimensional Phylogenetic Likelihoods." _Mol. Biol. Evol._ 35(1):242-246. https://doi.org/10.1093/molbev/msx253 [↩](#cite-7)
