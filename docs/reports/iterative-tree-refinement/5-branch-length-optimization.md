# Chapter 5: Per-edge branch length optimization

[Back to index](_index.md) | Previous: [Chapter 4: Ancestral reconstruction](4-ancestral-reconstruction.md) | Next: [Chapter 6: Zero-length branches](6-zero-length-branches.md)

## The inner optimization problem

Given a phylogenetic tree with fixed topology and fixed ancestral state distributions at every node (from [Chapter 4](4-ancestral-reconstruction.md)), find the branch length `t` for a single edge that maximizes the likelihood. This is a 1D scalar optimization problem -- the innermost layer of the tree refinement architecture.

The per-site conditional likelihood for one edge, using the eigendecomposition from [Chapter 2](2-substitution-models.md), is:

```
L_i(t) = sum_c k_{ic} * exp(lambda_c * t)
```

where `k_{ic}` are precomputed coefficients (depending on endpoint messages, not on t) and `lambda_c` are the model eigenvalues. The total log-likelihood is:

```
log L(t) = sum_i log L_i(t)
```

The optimizer maximizes this 1D function.

## Newton-Raphson (v1)

v1 uses Newton-Raphson with analytical derivatives from the eigendecomposition:

```
t_new = t - f'(t) / f''(t)
```

The derivatives are:

```
f'(t)  = sum_i [sum_c k_{ic} * lambda_c * exp(lambda_c * t)] / L_i(t)
f''(t) = sum_i {[sum_c k_{ic} * lambda_c^2 * exp(lambda_c * t)] / L_i(t) - [f'_i(t)]^2}
```

No numerical differentiation needed. The coefficients `k_{ic}` are computed once per edge (from `msg_to_parent` and `msg_to_child`) and reused across all Newton iterations.

Newton-Raphson converges quadratically near the optimum (each iteration doubles the number of correct digits). The step is clamped to `[-1.0, current_bl]` to prevent negative branch lengths.

**When Newton fails:** If the second derivative `f''(t) >= 0`, the surface is not concave and Newton would move toward a minimum. v1 falls back to a **grid search**: evaluate the likelihood on 100 logarithmically spaced branch lengths (`geomspace`) from `0.1 / L` to `1.5 * current_bl + 1/L`, and pick the maximum.

v1 code: `run_optimize_mixed()` in [`packages/treetime/src/commands/optimize/optimize_unified.rs#L533`](../../../packages/treetime/src/commands/optimize/optimize_unified.rs#L533). Up to 10 inner iterations per edge.

## Brent's method (v0)

v0 uses Brent's method (Brent 1973) -- a derivative-free bracket-based optimizer combining bisection, secant, and inverse quadratic interpolation. Brent's method guarantees convergence within a bracket `[a, b]` without requiring derivatives.

v0 applies Brent in **sqrt(t) space**: the variable `u = sqrt(t)` spreads out the near-zero region (where the likelihood surface is steep) and compresses the large-t region (where it is flat). This improves numerical stability for short branches.

v0 code: `optimal_marginal_branch_length()` calls `gtr.optimal_t_compressed()` which uses `scipy.optimize.minimize_scalar(method='brent')` at [`packages/legacy/treetime/treetime/gtr.py#L816-L920`](../../../packages/legacy/treetime/treetime/gtr.py#L816-L920). v0 also has an `fminbound` fallback at lines 904-910 for cases where Brent's method does not converge.

## Unimodality and local optima

For JC69 and F81 (all binary symmetric models), the per-edge log-likelihood has at most one stationary point -- it is **unimodal**. Any optimization method (Newton, Brent, grid search) converges to the global optimum (Dinh and Matsen 2017).

For K2P (Kimura 2-parameter) and more complex models (HKY, GTR), the likelihood **can have multiple local maxima**. Dinh and Matsen (2017) constructed explicit counterexamples for K2P and proved that the space of rescaled 1D likelihood functions is dense in all continuous non-negative functions on `[0, infinity)`. Any shape is possible.

Practical implication: for GTR models, Newton-Raphson from a single starting point can converge to a local optimum. The grid search fallback mitigates this by sampling the surface broadly. The initial guess ([Chapter 8](8-initial-estimation.md)) provides a starting point near the global optimum.

## Initial branch length estimation

Before the optimization loop starts, each branch needs an initial length estimate. The quality of this estimate affects convergence speed and correctness under multimodal likelihoods.

### v0: Hamming distance with Brent bracket

v0 uses the Hamming distance (count of positions differing between parent and child MAP sequences) divided by the alignment length as the initial estimate:

```
initial_bl = max(one_mutation, hamming_distance / alignment_length)
```

This seeds the bracket for Brent's method. Since Brent is bracket-based, it does not need a precise starting point -- it searches within `[0, 4.0]`.

### v1: substitution count from marginal profiles

v1 counts substitutions from the marginal reconstruction and divides by the effective alignment length (positions where both parent and child have canonical, non-gap states):

```
initial_bl = #subs / effective_length
```

The substitution count comes from `edge_subs()` -- MAP state comparison for dense partitions, marginal-reconstructed states for sparse partitions. Gap positions are excluded from both numerator and denominator.

v1 code: `initial_guess_mixed()` in [`packages/treetime/src/commands/optimize/optimize_unified.rs#L729`](../../../packages/treetime/src/commands/optimize/optimize_unified.rs#L729).

## Tool comparison

| Aspect               | v0                         | v1                    | RAxML-NG             | IQ-TREE              | PhyML              |
| -------------------- | -------------------------- | --------------------- | -------------------- | -------------------- | ------------------ |
| Per-edge method      | Brent in sqrt(t)           | Newton-Raphson        | Newton-Raphson       | NR + bisection       | Cubic spline on f' |
| Derivatives          | None (function only)       | 1st + 2nd analytical  | 1st + 2nd analytical | 1st + 2nd analytical | 1st only           |
| Non-concave fallback | N/A (bracket-based)        | Grid search (100 pts) | Step = -f/\|df\|     | Bisection step       | Bracket shrinking  |
| Max inner iterations | scipy default              | 10                    | 30                   | 100                  | 1020               |
| Branch bounds        | [0, 4.0]                   | [0, inf)              | [1e-6, 100]          | [1e-6, 10]           | [1e-8, 100]        |
| Near-zero handling   | sqrt(t) reparameterization | Derivative-sign check | z=exp(-t) reparam    | Clamp to min         | Geometric bracket  |
| Initial guess        | Hamming / L                | #subs / effective_L   | Tree builder value   | Tree builder value   | Tree builder value |

The key v0/v1 difference: v0 uses a derivative-free method (Brent) that handles non-concavity via bracketing but requires more function evaluations. v1 uses Newton-Raphson which converges faster when the surface is concave, with grid search fallback for non-concave cases.

RAxML, IQ-TREE, and PhyML do not compute initial guesses from sequence divergence -- they use the tree builder's branch lengths as starting points. Their reliability comes from per-branch safeguards (step clamping, bisection fallback, hard bounds), not from initial guess quality.

## References

- Felsenstein, J. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _J. Mol. Evol._ 17:368-376. https://doi.org/10.1007/BF01734359
- Brent, R. P. 1973. _Algorithms for Minimization without Derivatives._ Prentice-Hall. https://maths-people.anu.edu.au/~brent/pd/rpb011i.pdf
- Dinh, V. C., and F. A. Matsen IV. 2017. "The Shape of the One-Dimensional Phylogenetic Likelihood Function." _Ann. Appl. Prob._ 27(3):1646-1677. https://doi.org/10.1214/16-AAP1240
- Nocedal, J., and S. J. Wright. 2006. _Numerical Optimization._ 2nd ed. Springer. https://doi.org/10.1007/978-0-387-40065-5
- Stamatakis, A. 2014. "RAxML Version 8." _Bioinformatics_ 30(9):1312-1313. https://doi.org/10.1093/bioinformatics/btu033
- Nguyen, L.-T., H. A. Schmidt, A. von Haeseler, and B. Q. Minh. 2015. "IQ-TREE." _Mol. Biol. Evol._ 32(1):268-274. https://doi.org/10.1093/molbev/msu300
- Guindon, S., and O. Gascuel. 2003. "A Simple, Fast, and Accurate Algorithm." _Syst. Biol._ 52(5):696-704. https://doi.org/10.1080/10635150390235520
