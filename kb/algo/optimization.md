# Numerical Optimization Algorithms

[Back to index](README.md)

## Newton-Raphson for Branch Length Optimization

Per-edge branch length optimization using Newton's method with analytical first and second derivatives of the log-likelihood. Three parameterizations trade off singularity handling near zero:

- t-space (`newton_inner`): update `t_new = t - clamp(f'/f'', -1.0, t)`. Indel Hessian $-k/t^2$ dominates on short branches.
- sqrt(t)-space (`newton_sqrt_inner`): chain rule with $s = \sqrt{t}$. Reduces indel singularity from $O(1/t^2)$ to $O(1/t)$.
- ln(t)-space (`newton_log_inner`): chain rule with $u = \ln(t)$. Eliminates indel singularity entirely (bounded curvature $-\mu t$). Correct step clamping: upper $u - u_{\min}$, lower $-\ln(1 + 1/t)$.

A non-negative initial second derivative triggers a 100-point log-spaced grid search. If the second derivative becomes non-negative after a Newton step, the iteration stops and returns its current candidate.

v1: [`packages/treetime/src/optimize/method_newton.rs`](../../packages/treetime/src/optimize/method_newton.rs).

v0 uses Brent's method (`scipy.optimize.minimize_scalar`) in $\sqrt{t}$ space with a Hamming-distance bracket. v1 ships six methods (Newton and Brent in $t$, $\sqrt{t}$, and $\ln(t)$ coordinates) selectable via `--opt-method`; `brent-sqrt` is the default. It uses the same square-root parameterization and corresponding substitution-likelihood objective form as v0, but the solver, bounds, tolerances, regularization, indel term, and resulting branch lengths are not exactly equivalent. The Newton variants compute analytical derivatives from cached eigenvalue-space coefficients. The implementation evaluates the second derivative in a centered eigenvalue form to reduce cancellation relative to subtracting two independently accumulated moments; its signed weights are not a posterior distribution or Welford recurrence. See [feature inventory](../features/optimize.md) for parity evidence and [intentional change](../decisions/optimize-newton-raphson-per-edge.md) for the per-method rationale.

The method is described in <a id="cite-1"></a>[Nocedal and Wright 2006](https://doi.org/10.1007/978-0-387-40065-5) [[1](#ref-1)], Chapter 2, and <a id="cite-2"></a>[Felsenstein 2003](https://doi.org/10.1007/978-0-387-21337-7) [[2](#ref-2)], Chapter 16.

---

## Brent's Method for Branch Length Optimization

Per-edge branch length optimization using Brent's derivative-free method (<a id="cite-3"></a>[Brent 1973](https://maths-people.anu.edu.au/brent/pub/pub011.html) [[3](#ref-3)]) via `argmin::BrentOpt`. Three coordinate transformations change the bounds, curvature, and tolerance interpretation:

- t-space (`brent_inner`): bracket $[\text{min\_bl}, \text{upper}]$, cost function evaluates $-\ell(t)$. Baseline derivative-free method.
- sqrt(t)-space (`brent_sqrt_inner`): bracket $[\sqrt{\text{min\_bl}}, \sqrt{\text{upper}}]$, cost function evaluates $-\ell(s^2)$. This is the default and uses v0's parameterization. Tolerance $\epsilon_s$ in $s$-space maps locally to $t$-space precision $\approx 2s^* \epsilon_s$.
- ln(t)-space (`brent_log_inner`): bracket $[\ln(\text{min\_bl}), \ln(\text{upper})]$, cost function evaluates $-\ell(e^u)$. Tolerance $\epsilon_u$ is locally a relative tolerance because $dt/t \approx du$.

v1: [`packages/treetime/src/optimize/method_brent.rs`](../../packages/treetime/src/optimize/method_brent.rs).

v0: `optimal_t_compressed()` at [`packages/legacy/treetime/treetime/gtr.py#L816-L920`](../../packages/legacy/treetime/treetime/gtr.py#L816-L920). Uses `scipy.optimize.minimize_scalar(method='brent')` in sqrt(t) space with bracket `[-sqrt(MAX_BRANCH_LENGTH), sqrt(hamming_distance), sqrt(MAX_BRANCH_LENGTH)]`.

Distinct from existing Brent entries for clock root optimization, coalescent Tc optimization, and polytomy resolution. Those are different optimization targets using the same `BrentOpt` solver.

See also <a id="cite-4"></a>[Press et al. 2007](https://doi.org/10.1017/CBO9780511811340) [[4](#ref-4)], Section 10.3.

---

## Eigendecomposition-Based Likelihood

Precomputes eigenvector coefficients for each edge, enabling efficient per-branch log-likelihood and derivative evaluation without repeated matrix exponential computation. For an eigendecomposition

$$Q=V\operatorname{diag}(\lambda)V^{-1},$$

the transition matrix is

$$P(t)=V\operatorname{diag}(e^{\lambda_c t})V^{-1}.$$

If the substitution rate has not already been absorbed into $Q$ or the branch coordinate, the exponent is $\mu\lambda_c t$. Profile-eigenvector products such as `msg.dot(V)` and `msg.dot(V_inv.T)` are independent of branch length and can be cached once per edge.

v1 dense: [`packages/treetime/src/optimize/dense_eval.rs`](../../packages/treetime/src/optimize/dense_eval.rs).
v1 sparse: [`packages/treetime/src/optimize/sparse_eval.rs`](../../packages/treetime/src/optimize/sparse_eval.rs).

Sparse evaluation handles variable positions individually and aggregates conserved positions by fixed-state counts.

The eigendecomposition approach follows <a id="cite-5"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[5](#ref-5)] and <a id="cite-6"></a>[Yang 2006](https://doi.org/10.1093/acprof:oso/9780198567028.001.0001) [[6](#ref-6)], Chapter 4.

---

## Poisson Indel Contribution

Adds a Poisson term to branch-length optimization. Here $k$ is the current number of net contiguous indel annotation records on one edge, $t$ is that edge's branch length, $\mu$ is the global record rate, and $e$ indexes tree edges. The per-edge contribution is $\ell(t) = k \ln(\mu t) - \mu t - \ln(k!)$. Derivatives $k/t - \mu$ and $-k/t^2$ enter the Newton step alongside substitution derivatives. Because annotation composition can merge or cancel records, $k$ is not generally a count of historical indel events. The estimate $\hat{\mu} = \sum_e k_e / \sum_e t_e$ is computed before the internal iterations of one `run_optimize_loop()` invocation and then held fixed for its tree-level objective and edge updates.

v1: [`packages/treetime/src/optimize/indel.rs`](../../packages/treetime/src/optimize/indel.rs).

v0: not implemented; v0 ignores indels in this branch-length likelihood.

This is a v1-only feature. See [indel models report](../reports/indel-models/README.md) for the full catalog of indel modeling approaches, [intentional change](../decisions/optimize-indel-contribution-to-likelihood.md), [design doc](../_raw/optimize.md), and [alternatives proposal](../proposals/optimize-indel-model-alternatives.md).

---

## Piecewise Linear Interpolation (Uniform Grid)

O(1) interval lookup via `floor((x - x_min) / dx)` for uniformly spaced grids. Used throughout the distribution system for evaluating discretized probability distributions and branch length likelihoods on fixed grids.

v1: [`packages/treetime-grid/src/grid_fn.rs#L279-L351`](../../packages/treetime-grid/src/grid_fn.rs#L279-L351).

Standard treatment in <a id="cite-7"></a>[Burden and Faires 2010](https://doi.org/10.1007/978-3-319-07671-3) [[7](#ref-7)], Chapter 3.

---

## Piecewise Linear Interpolation (Non-Uniform Grid)

O(log n) binary search interval lookup for non-uniformly spaced grids. Used for the skyline coalescent Tc(t) function where grid points are placed at coalescent event times rather than on a uniform grid.

v1: [`packages/treetime-grid/src/interp_nonuniform.rs#L25-L56`](../../packages/treetime-grid/src/interp_nonuniform.rs#L25-L56).

---

## Exponential Damping for Outer-Loop Convergence

Blends optimized branch lengths with previous values using iteration-dependent weights: `bl = bl_new * (1 - damping^(i+1)) + bl_old * damping^(i+1)`. Early iterations take conservative steps; later iterations approach the full update. This is intended to moderate oscillation in the alternating marginal-reconstruction and branch-length-update cycle.

v1: [`packages/treetime/src/commands/optimize/run.rs`](../../packages/treetime/src/commands/optimize/run.rs) `apply_damping()`.

v0: `optimize_tree_marginal()` at [`packages/legacy/treetime/treetime/treeanc.py#L1297-L1360`](../../packages/legacy/treetime/treetime/treeanc.py#L1297-L1360) with `damping=0.75` default.

Based on <a id="cite-8"></a>[Sagulenko, Puller, and Neher 2018](https://doi.org/10.1093/ve/vex042) [[8](#ref-8)].

---

## Three-Condition Convergence Check

The optimize outer loop uses three stopping conditions in `run_optimize_loop()` ([`packages/treetime/src/optimize/run_loop.rs`](../../packages/treetime/src/optimize/run_loop.rs)), where $L_i$ is the joint log-likelihood at iteration $i$ and $\mathit{dp}$ is the convergence threshold. In indel-bearing runs, $L_i$ includes both substitution likelihood from `update_marginal()` and the tree-level Poisson contribution evaluated with the fixed rate for that loop invocation:

Aggregate likelihood results and histories use `LogLh`; per-site coefficients, derivatives, branch lengths, and solver coordinates remain `f64`. Subtraction between two `LogLh` values produces the scalar convergence delta used below. See [kb/decisions/typed-log-likelihood-values.md](../decisions/typed-log-likelihood-values.md).

- Small change: $|L_i - L_{i-1}| < \mathit{dp}$. Because this uses an absolute difference, it can also accept a small deterioration.
- Two-cycle: $|L_i - L_{i-2}| < \mathit{dp}$ for $i\ge2$.
- Decline after a best iteration: $i\ge2$, $L_i<L_{i-1}$, and $L_{i-1}\ge L_{\mathrm{best}}$; the saved branch lengths are restored and iteration stops.

The damping schedule applies a floor: $d_i = \max(d^{i+1},\, \texttt{DAMPING\_FLOOR})$ where $\texttt{DAMPING\_FLOOR} = 0.01$, preventing fully undamped late iterations that amplify the sparse discrete jump. A NaN/Inf guard breaks immediately on numerical instability.

v1 defaults: `max_iter=10`, `dp=0.1`, matching v0's `optimize_tree_marginal()`.

v0: uses a signed check (`deltaLH < LHtol`) which conflates convergence with worsening. See [v0 erratum](../v0-errata/optimize-signed-convergence-check.md).

Two-cycles have been observed in sparse optimization. The implemented loop combines coordinate-wise edge optimization, damping, topology edits, profile updates, and a fixed indel-rate estimate, so it has not been established as an exact EM algorithm with a monotonicity guarantee.

IQ-TREE's per-round monotonicity check described in <a id="cite-9"></a>[Minh et al. 2020](https://doi.org/10.1093/molbev/msaa015) [[9](#ref-9)].

---

## Edge Collapse (Topology Cleanup)

Canonical edge-collapse operation for zero-length or near-zero-length internal edges. Reparents children to the grandparent, updates branch lengths (summing parent and child), and cleans up stale partition entries for both sparse (sub composition) and dense representations. Shared across the optimize and prune commands.

v1: `collapse_edge()` in [`packages/treetime/src/optimize/topology/collapse.rs`](../../packages/treetime/src/optimize/topology/collapse.rs).

v0: inline in `prune_short_branches()` at [`packages/legacy/treetime/treetime/treeanc.py#L1475-L1496`](../../packages/legacy/treetime/treetime/treeanc.py#L1475-L1496).

---

## Forward-Pass Zero-Divisor Clamping

The sparse and dense forward-pass marginal divisions floor zero divisors at `f64::MIN_POSITIVE`. This prevents NaN propagation, but also turns a structural zero into a finite contribution and therefore changes impossible-state semantics. `normalize_inplace` returns a uniform distribution for zero-sum or non-finite rows. The scientific policy for these degenerate cases remains unresolved.

v1 sparse: [`packages/treetime/src/partition/marginal_passes.rs`](../../packages/treetime/src/partition/marginal_passes.rs).
v1 dense: [`packages/treetime/src/partition/marginal_dense.rs`](../../packages/treetime/src/partition/marginal_dense.rs).

v0: no explicit guard; relies on NumPy's inf/nan propagation behavior.

---

## Zero-Length Branch Detection (Derivative Sign)

For contributions whose model is flagged as unimodal, evaluates the sign of $d\log L/dt$ at zero. With $L_i(t)=\sum_c k_{ic}\exp(\lambda_ct)$, the per-site derivative at zero is $(\sum_c k_{ic}\lambda_c)/(\sum_c k_{ic})$. A negative total derivative implies that zero is optimal only under the unimodality restriction; richer time-reversible models can have multimodal one-dimensional likelihoods.

v1: [`packages/treetime/src/optimize/zero_boundary.rs`](../../packages/treetime/src/optimize/zero_boundary.rs) `is_zero_branch_optimal()`.

v0 does not have an equivalent analytical check. v0 uses `prune_short_branches()` with a compound threshold-and-probability criterion instead.

The model restriction is consistent with the one-dimensional likelihood analysis of <a id="cite-10"></a>[Dinh and Matsen 2017](https://doi.org/10.1214/16-AAP1240) [[10](#ref-10)].

---

## Zero-Length Branch Pruning

Collapses internal edges whose optimal length is zero or near-zero, reparenting children to the grandparent. Creates polytomies. v0 uses a compound criterion: `branch_length < 0.1 * one_mutation AND prob_t(parent_seq, child_seq, 0) > 0.1`.

v0: `prune_short_branches()` at [`packages/legacy/treetime/treetime/treeanc.py#L1475-L1496`](../../packages/legacy/treetime/treetime/treeanc.py#L1475-L1496).

v1 collapses internal edges driven to exactly zero during optimization. It does not implement v0's post-convergence near-zero threshold plus probability criterion. The prune command separately provides `--prune-short` and `--prune-empty`.

---

## Shared-Mutation Polytomy Merging

Scans children of polytomy nodes for shared substitutions. When two siblings carry identical mutations, they are grouped under a new internal node whose branch length is the Jukes-Cantor 1969 correction of the pooled p-distance `#shared_mutations / alignment_length` (see [`jukes_cantor_distance()`](../../packages/treetime/src/gtr/jc_distance.rs)). Shared mutations move to the new parent edge; remaining unique mutations stay on child edges. See [the corresponding intentional change](../decisions/prune-merge-jukes-cantor-branch-length.md) for the rationale for correcting the raw ratio specified in `../_raw/optimize.md`.

v1: [`packages/treetime/src/optimize/topology/merge_shared_mutations.rs`](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs) `merge_shared_mutation_branches()`. Invoked by `prune --merge-shared-mutations` and by the `optimize` topology-cleanup pre-step before each per-edge optimization round.

v0: No formal implementation. Design doc describes "ad-hoc scripts" in nextstrain pathogen pipelines.

---

## Greedy Temporal Polytomy Resolution

For each polytomy, computes pairwise likelihood gain from merging children under a new intermediate node. Uses Brent optimization over the time domain with `zero_branch_slope = mu * L` penalty for the new zero-mutation branch. Greedily picks the best pair above `resolution_threshold` (0.05). O(n^2) per polytomy.

v1: [`packages/treetime/src/timetree/optimization/polytomy.rs`](../../packages/treetime/src/timetree/optimization/polytomy.rs) `resolve_polytomies()`.

v0: `_poly()` at [`packages/legacy/treetime/treetime/treetime.py#L713-L870`](../../packages/legacy/treetime/treetime/treetime.py#L713-L870). v0 distinguishes "stretched" (`mutation_length < clock_length`) from "compressed" children and by default only resolves stretched ones.

---

## Stochastic Coalescent Polytomy Resolution

Simulates a backward-in-time coalescent process. Branches without mutations are "ready to coalesce"; branches with mutations must have mutations removed stochastically. Randomly pairs ready branches for merging. Recommended for large polytomies where greedy mode produces caterpillar-like subtrees.

v0: `generate_subtree()` at [`packages/legacy/treetime/treetime/treetime.py#L872-L1010`](../../packages/legacy/treetime/treetime/treetime.py#L872-L1010).

v1: Not implemented. Tracked: `N-timetree-stochastic-polytomy-unimplemented.md`.

---

## References

- <a id="ref-1"></a>Nocedal, Jorge, and Stephen J. Wright. 2006. _Numerical Optimization._ 2nd ed. Springer. ISBN 978-0-387-30303-1. https://doi.org/10.1007/978-0-387-40065-5 [↩](#cite-1)
- <a id="ref-2"></a>Felsenstein, Joseph. 2003. _Inferring Phylogenies._ Sinauer Associates. ISBN 978-0-87893-177-4. [↩](#cite-2)
- <a id="ref-3"></a>Brent, Richard P. 1973. _Algorithms for Minimization Without Derivatives._ Prentice-Hall. ISBN 0-13-022335-2. https://maths-people.anu.edu.au/brent/pub/pub011.html [↩](#cite-3)
- <a id="ref-4"></a>Press, William H., Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery. 2007. _Numerical Recipes: The Art of Scientific Computing._ 3rd ed. Cambridge University Press. ISBN 978-0-521-88068-8. https://doi.org/10.1017/CBO9780511811340 [↩](#cite-4)
- <a id="ref-5"></a>Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _Journal of Molecular Evolution_ 17(6):368-376. https://doi.org/10.1007/BF01734359 [↩](#cite-5)
- <a id="ref-6"></a>Yang, Ziheng. 2006. _Computational Molecular Evolution._ Oxford University Press. ISBN 978-0-19-856702-8. [↩](#cite-6)
- <a id="ref-7"></a>Burden, Richard L., and J. Douglas Faires. 2010. _Numerical Analysis._ 9th ed. Brooks/Cole. ISBN 978-0-538-73351-9. [↩](#cite-7)
- <a id="ref-8"></a>Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4(1):vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-8)
- <a id="ref-9"></a>Minh, Bui Quang, Heiko A. Schmidt, Olga Chernomor, Dominik Schrempf, Michael D. Woodhams, Arndt von Haeseler, and Robert Lanfear. 2020. "IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era." _Molecular Biology and Evolution_ 37(5):1530-1534. https://doi.org/10.1093/molbev/msaa015 [↩](#cite-9)
- <a id="ref-10"></a>Dinh, Vu, and Frederick A. Matsen IV. 2017. "The Shape of the One-Dimensional Phylogenetic Likelihood Function." _The Annals of Applied Probability_ 27(1):1-17. https://doi.org/10.1214/16-AAP1240 [↩](#cite-10)

---

## File Index

| File                                                                                                                               | Algorithms                                                                             |
| ---------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------- |
| [`packages/treetime/src/optimize/`](../../packages/treetime/src/optimize/)                                       | Newton-Raphson, Brent, grid search, likelihood eval, damping, convergence, zero-detect |
| [`packages/treetime/src/prune/`](../../packages/treetime/src/prune/)                                             | Shared-mutation merging                                                                |
| [`packages/treetime/src/timetree/optimization/`](../../packages/treetime/src/timetree/optimization/)             | Greedy temporal polytomy resolution                                                    |
| [`packages/treetime/src/optimize/topology/`](../../packages/treetime/src/optimize/topology/) | Edge collapse (shared)                                                                 |
| [`packages/treetime/src/partition/`](../../packages/treetime/src/partition/)                         | Forward-pass zero-divisor clamping, normalize_inplace                                  |
| [`packages/treetime-grid/src/`](../../packages/treetime-grid/src/)                                                                 | Interpolation (uniform, non-uniform)                                                   |
