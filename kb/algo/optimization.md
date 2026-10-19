# Numerical Optimization Algorithms

## Newton-Raphson for Branch Length Optimization

Per-edge branch length optimization using Newton's method with analytical first and second derivatives of the log-likelihood, in sqrt(t) space (`solve_newton`): chain rule with $s = \sqrt{t}$, update $s_{\text{new}} = s - \operatorname{clamp}(\ell'_s/\ell''_s,\ s - \sqrt{s^2 + 1},\ s)$. The lower step bound limits one step to a $t$ increase of 1. The parameterization reduces the Poisson indel Hessian singularity from $O(1/t^2)$ in $t$ space to $O(1/t)$.

The search keeps a bracket $[s_{lo}, s_{hi}]$ in $s$-space that the sign of the slope narrows after each evaluation. A Newton step is taken when the curvature is negative and the step stays inside the bracket; otherwise the iteration bisects the bracket away from the current point. The iteration reports convergence only at a concave maximum within tolerance or when the bracket is narrower than the tolerance, and returns the best evaluated point. After 50 iterations it returns the best point with the status iteration limit, which raises the `branch-lengths-not-converged` run warning.

Newton's method converges only from a start inside its basin of attraction. Input branch lengths usually satisfy this. Edges with missing or invalid lengths, or all edges with `--branch-length-initial-guess=always`, start from the observed substitution count divided by the edge's effective (non-gap) length, or from the indel count over the indel rate on an edge with indels only. v1: `fn guess_initial_branch_lengths()` in [`packages/treetime/src/optimize/dispatch.rs#L167`](../../packages/treetime/src/optimize/dispatch.rs#L167).

v1: [`packages/treetime/src/optimize/method_newton.rs`](../../packages/treetime/src/optimize/method_newton.rs).

v0 uses Brent's method (`scipy.optimize.minimize_scalar`) in $\sqrt{t}$ space with a Hamming-distance bracket. v1 ships two methods (Brent and Newton in $\sqrt{t}$ coordinates) selectable via `--opt-method`; `brent` is the default. It uses the same square-root parameterization and corresponding substitution-likelihood objective form as v0, but the solver, bounds, tolerances, regularization, indel term, and resulting branch lengths are not exactly equivalent. The Newton method computes analytical derivatives from cached eigenvalue-space coefficients. The implementation evaluates the second derivative in a centered eigenvalue form to reduce cancellation relative to subtracting two independently accumulated moments; its signed weights are not a posterior distribution or Welford recurrence. See [feature inventory](../features/optimize.md) for parity evidence and [intentional change](../decisions/optimize-newton-raphson-per-edge.md) for the per-method rationale.

The method is described in <a id="cite-1"></a>[Nocedal and Wright 2006](https://doi.org/10.1007/978-0-387-40065-5) [[1](#ref-1)], Chapter 2, and <a id="cite-2"></a>[Felsenstein 2003](https://www.oup.com.au/books/higher-education/biology/9780878931774) [[2](#ref-2)], Chapter 16.

---

## Brent's Method for Branch Length Optimization

Per-edge branch length optimization using Brent's derivative-free method (<a id="cite-3"></a>[Brent 1973](https://maths-people.anu.edu.au/brent/pub/pub011.html) [[3](#ref-3)]) via `argmin::BrentOpt`, in sqrt(t) space (`solve_brent`): bracket $[\sqrt{\text{min\_bl}}, \sqrt{\text{upper}}]$, cost function evaluates $-\ell(s^2)$. This is the default and uses v0's parameterization. Tolerance $\epsilon_s$ in $s$-space maps locally to $t$-space precision $\approx 2s^* \epsilon_s$.

v1: [`packages/treetime/src/optimize/method_brent.rs`](../../packages/treetime/src/optimize/method_brent.rs).

v0: `optimal_t_compressed()` at [`packages/legacy/treetime/treetime/gtr.py#L816-L920`](../../packages/legacy/treetime/treetime/gtr.py#L816-L920). Uses `scipy.optimize.minimize_scalar(method='brent')` in sqrt(t) space with bracket `[-sqrt(MAX_BRANCH_LENGTH), sqrt(hamming_distance), sqrt(MAX_BRANCH_LENGTH)]`.

Distinct from existing Brent entries for clock root optimization and coalescent Tc optimization. Those are different optimization targets using the same `BrentOpt` solver.

See also <a id="cite-4"></a>[Press et al. 2007](https://doi.org/10.1017/CBO9780511811340) [[4](#ref-4)], Section 10.3.

---

## Eigendecomposition-Based Likelihood

Precomputes eigenvector coefficients for each edge, enabling efficient per-branch log-likelihood and derivative evaluation without repeated matrix exponential computation. For an eigendecomposition

$$Q=V\operatorname{diag}(\lambda)V^{-1},$$

the transition matrix is

$$P(t)=V\operatorname{diag}(e^{\mu\lambda_c t})V^{-1}.$$

$Q$ here is the normalized rate matrix of the GTR, and $\mu$ is its overall rate. The objective evaluates the exponent $\mu\lambda_c t$ with the eigenvalues of $\mu Q$ (`GTR::rate_eigvals()`), the same transition matrix that the marginal passes use, so a change of one branch length changes the objective and the tree log likelihood by the same amount. Profile-eigenvector products such as `msg.dot(V)` and `msg.dot(V_inv.T)` are independent of branch length and can be cached once per edge.

v1 dense: [`packages/treetime/src/optimize/dense_eval.rs`](../../packages/treetime/src/optimize/dense_eval.rs).
v1 sparse: [`packages/treetime/src/optimize/sparse_eval.rs`](../../packages/treetime/src/optimize/sparse_eval.rs).

Sparse evaluation handles variable positions individually and aggregates conserved positions by fixed-state counts.

The eigendecomposition approach follows <a id="cite-5"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[5](#ref-5)] and <a id="cite-6"></a>[Yang 2006](https://doi.org/10.1093/acprof:oso/9780198567028.001.0001) [[6](#ref-6)], Chapter 4.

---

## Poisson Indel Contribution

Adds a Poisson term to branch-length optimization. Here $k$ is the current number of net contiguous indel annotation records on one edge, $t$ is that edge's branch length, $\mu$ is the global record rate, and $e$ indexes tree edges. The per-edge contribution is $\ell(t) = k \ln(\mu t) - \mu t - \ln(k!)$. Derivatives $k/t - \mu$ and $-k/t^2$ enter the Newton step alongside substitution derivatives. Annotation composition can merge or cancel records, so $k$ counts current records instead of historical indel events. The estimate $\hat{\mu} = \sum_e k_e / \sum_e t_e$ is computed before the internal iterations of one `run_optimize_loop()` invocation and then held fixed for its tree-level objective and edge updates.

v1: [`packages/treetime/src/optimize/indel.rs`](../../packages/treetime/src/optimize/indel.rs).

v0: not implemented; v0 ignores indels in this branch-length likelihood.

This is a v1-only feature. See the indel models report ([kb/reports/indel-models/1-introduction.md](../reports/indel-models/1-introduction.md)) for the full catalog of indel modeling approaches, [intentional change](../decisions/optimize-indel-contribution-to-likelihood.md) and [alternatives proposal](../proposals/optimize-indel-model-alternatives.md).

---

## Piecewise Linear Interpolation (Uniform Grid)

O(1) interval lookup via `floor((x - x_min) / dx)` for uniformly spaced grids. Used throughout the distribution system for evaluating discretized probability distributions and branch length likelihoods on fixed grids.

v1: [`packages/treetime-grid/src/grid_fn.rs#L274-L345`](../../packages/treetime-grid/src/grid_fn.rs#L274-L345).

Standard treatment in <a id="cite-7"></a>[Burden and Faires 2010](https://doi.org/10.1007/978-3-319-07671-3) [[7](#ref-7)], Chapter 3.

---

## Piecewise Linear Interpolation (Non-Uniform Grid)

O(log n) binary search interval lookup for non-uniformly spaced grids. Used for the skyline coalescent Tc(t) function where grid points are placed at coalescent event times rather than on a uniform grid.

v1: [`packages/treetime-grid/src/interp_nonuniform.rs#L25-L56`](../../packages/treetime-grid/src/interp_nonuniform.rs#L25-L56).

---

## Exponential Damping for Outer-Loop Convergence

Blends optimized branch lengths with previous values using iteration-dependent weights: `bl = bl_new * (1 - damping^(i+1)) + bl_old * damping^(i+1)`. Early iterations take conservative steps; later iterations approach the full update. This is intended to moderate oscillation in the alternating marginal-reconstruction and branch-length-update cycle.

v1: `fn apply_damping()` in [`packages/treetime/src/optimize/iteration.rs#L7`](../../packages/treetime/src/optimize/iteration.rs#L7).

v0: `optimize_tree_marginal()` at [`packages/legacy/treetime/treetime/treeanc.py#L1297-L1360`](../../packages/legacy/treetime/treetime/treeanc.py#L1297-L1360) with `damping=0.75` default.

Based on <a id="cite-8"></a>[Sagulenko, Puller, and Neher 2018](https://doi.org/10.1093/ve/vex042) [[8](#ref-8)].

---

## Three-Condition Convergence Check

The optimize outer loop uses three stopping conditions in `run_optimize_loop()` ([`packages/treetime/src/optimize/run_loop.rs`](../../packages/treetime/src/optimize/run_loop.rs)), where $L_i$ is the joint log-likelihood at iteration $i$ and $\mathit{dp}$ is the convergence threshold. In indel-bearing runs, $L_i$ includes both substitution likelihood from `update_marginal()` and the tree-level Poisson contribution evaluated with the fixed rate for that loop invocation:

Aggregate likelihood results and histories use `LogLh`; per-site coefficients, derivatives, branch lengths, and solver coordinates remain `f64`. Subtraction between two `LogLh` values produces the scalar convergence delta used below.

- Small change: $|L_i - L_{i-1}| < \mathit{dp}$. Because this uses an absolute difference, it can also accept a small deterioration.
- Two-cycle: $|L_i - L_{i-2}| < \mathit{dp}$ for $i\ge2$.
- Decline after a best iteration: $i\ge2$, $L_i<L_{i-1}$, and $L_{i-1}\ge L_{\mathrm{best}}$.

When the loop exhausts `max_iter`, it evaluates the branch lengths of the last update once more. Whatever the stop reason, it returns the evaluated branch lengths with the highest $L_i$ since the last topology change. Simultaneous per-edge updates need not increase the likelihood: edges whose likelihood depends almost only on their sum, such as a path through nodes whose other subtrees carry no data, overshoot together, and $L_i$ can decline from the first iteration and never trigger the decline condition. Edges joined by a node with exactly two edges are optimized as one sum and do not overshoot ([kb/decisions/optimize-series-edges-joint-sum.md](../decisions/optimize-series-edges-joint-sum.md)).

With `--damping` above 0, the damping schedule applies a floor: $d_i = \max(d^{i+1},\, \texttt{DAMPING\_FLOOR})$ where $\texttt{DAMPING\_FLOOR} = 0.01$, preventing fully undamped late iterations that amplify the sparse discrete jump. The default `--damping` of 0 applies no damping ([kb/decisions/optimize-damping-disabled-by-default.md](../decisions/optimize-damping-disabled-by-default.md)). A NaN/Inf guard breaks immediately on numerical instability.

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

v1: `fn divide_out_value()` in [`packages/treetime/src/partition/marginal/shared/normalize.rs#L51-L53`](../../packages/treetime/src/partition/marginal/shared/normalize.rs#L51-L53), applied by the dense forward pass through `fn divide_out()` ([`packages/treetime/src/partition/marginal/shared/pass.rs`](../../packages/treetime/src/partition/marginal/shared/pass.rs)) and by the sparse forward pass through `fn divide_out_site()` ([`packages/treetime/src/partition/marginal/sparse/forward.rs`](../../packages/treetime/src/partition/marginal/sparse/forward.rs)).

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

Scans children of polytomy nodes for shared substitutions. When two siblings carry identical mutations, they are grouped under a new internal node whose branch length is the Jukes-Cantor 1969 correction of the pooled p-distance `#shared_mutations / alignment_length` (see [`jukes_cantor_distance()`](../../packages/treetime/src/gtr/jc_distance.rs)). Shared mutations move to the new parent edge; remaining unique mutations stay on child edges. See [the corresponding intentional change](../decisions/prune-merge-jukes-cantor-branch-length.md) for the rationale for correcting the raw ratio.

v1: [`packages/treetime/src/optimize/topology/merge_shared_mutations.rs`](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs) `merge_shared_mutation_branches()`. Invoked by `prune --merge-shared-mutations` and by the `optimize` topology-cleanup pre-step before each per-edge optimization round.

v0: No formal implementation. Design doc describes "ad-hoc scripts" in nextstrain pathogen pipelines.

---

## Stochastic Coalescent Polytomy Resolution

Simulates a backward-in-time coalescent process. Branches without remaining substitutions are ready to coalesce. Other branches must receive mutation events before they can coalesce. The event hazard is integrated across lineage arrivals and changes in the piecewise-constant coalescent rate.

v0: `generate_subtree()` at [`packages/legacy/treetime/treetime/treetime.py#L872-L1010`](../../packages/legacy/treetime/treetime/treetime.py#L872-L1010).

v1: [`packages/treetime/src/timetree/optimization/polytomy/`](../../packages/treetime/src/timetree/optimization/polytomy/). v1 uses consistent mutation weights and stops at the parent bound. See [kb/algo/timetree.md](timetree.md#polytomy-resolution) for the rates and state transitions.

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

| File                                                                                                 | Algorithms                                                                             |
| ---------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------- |
| [`packages/treetime/src/optimize/`](../../packages/treetime/src/optimize/)                           | Newton-Raphson, Brent, grid search, likelihood eval, damping, convergence, zero-detect |
| [`packages/treetime/src/prune/`](../../packages/treetime/src/prune/)                                 | Shared-mutation merging                                                                |
| [`packages/treetime/src/timetree/optimization/`](../../packages/treetime/src/timetree/optimization/) | Stochastic temporal polytomy resolution                                                |
| [`packages/treetime/src/optimize/topology/`](../../packages/treetime/src/optimize/topology/)         | Edge collapse (shared)                                                                 |
| [`packages/treetime/src/partition/`](../../packages/treetime/src/partition/)                         | Forward-pass zero-divisor clamping, normalize_inplace                                  |
| [`packages/treetime-grid/src/`](../../packages/treetime-grid/src/)                                   | Interpolation (uniform, non-uniform)                                                   |
