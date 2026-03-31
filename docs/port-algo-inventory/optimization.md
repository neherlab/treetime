# Numerical Optimization Algorithms

[Back to index](_index.md)

## Newton-Raphson for Branch Length Optimization

Per-edge branch length optimization using Newton's method with analytical first and second derivatives of the log-likelihood. The update rule is `t_new = t - clamp(f'/f'', -1.0, t)`, clamping the step to prevent negative branch lengths. When the second derivative is non-negative (likelihood surface is convex at the current point), the method falls back to a 100-point linear grid search over the branch length domain.

v1: [`packages/treetime/src/commands/optimize/optimize_unified.rs#L249-L268`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L249-L268).

v0 uses Brent's method (`scipy.optimize.minimize_scalar`) in sqrt(t) space with Hamming distance bracket instead of Newton's method. v1's analytical derivatives avoid the derivative-free overhead of Brent but require correct second-derivative computation. See [feature inventory](../port-feature-inventory/_index.md#7-branch-length-optimization) for parity details.

References:

- Nocedal & Wright. "Numerical Optimization." Chapter 2.
- Felsenstein (2003). "Inferring Phylogenies." Chapter 16 (branch length optimization in ML phylogenetics).

---

## Eigendecomposition-Based Likelihood

Precomputes eigenvector coefficients for each edge, enabling efficient per-branch log-likelihood and derivative evaluation without repeated matrix exponential computation. For a GTR model with rate matrix Q = V _ diag(lambda) _ V^-1, the transition probabilities at branch length t factor as `P(t) = V * diag(exp(lambda_i * t)) * V^-1`. The key insight is that profile-eigenvector dot products (`msg.dot(V)` and `msg.dot(V_inv.T)`) are branch-length-independent and can be cached once per edge.

v1 dense: [`packages/treetime/src/commands/optimize/optimize_dense_eval.rs`](../../packages/treetime/src/commands/optimize/optimize_dense_eval.rs).
v1 sparse: [`packages/treetime/src/commands/optimize/optimize_sparse_eval.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse_eval.rs).

The sparse path weights each site contribution by its multiplicity (number of identical columns in the alignment sharing that substitution pattern), reducing computation for conserved sequences.

References:

- Felsenstein (1981). "Evolutionary trees from DNA sequences." J Mol Evol, 17(6):368-376. doi:10.1007/BF01734359
- Yang (2006). "Computational Molecular Evolution." Chapter 4.

---

## Poisson Indel Contribution

Adds a Poisson indel log-likelihood term to per-edge branch length optimization. For $k$ observed indel events on a branch of length $t$ with global rate $\mu$: $\ell(t) = k \ln(\mu t) - \mu t - \ln(k!)$. Derivatives $k/t - \mu$ and $-k/t^2$ enter the Newton step alongside substitution derivatives. The rate $\hat{\mu} = \sum_e k_e / \sum_e t_e$ is estimated from the tree at each optimization round.

v1: [`packages/treetime/src/commands/optimize/optimize_indel.rs`](../../packages/treetime/src/commands/optimize/optimize_indel.rs).

v0: not implemented. v0 ignores indels in the likelihood, same as RAxML, IQ-TREE, PhyML.

This is a v1-only feature. See [intentional change](../port-intentional-changes/optimize-indel-contribution-to-likelihood.md) and [design doc](../algorithms/optimize.md).

---

## Piecewise Linear Interpolation (Uniform Grid)

O(1) interval lookup via `floor((x - x_min) / dx)` for uniformly spaced grids. Used throughout the distribution system for evaluating discretized probability distributions and branch length likelihoods on fixed grids.

v1: [`packages/treetime-grid/src/grid_fn.rs#L279-L351`](../../packages/treetime-grid/src/grid_fn.rs#L279-L351).

Reference: Burden & Faires. "Numerical Analysis." Chapter 3.

---

## Piecewise Linear Interpolation (Non-Uniform Grid)

O(log n) binary search interval lookup for non-uniformly spaced grids. Used for the skyline coalescent Tc(t) function where grid points are placed at coalescent event times rather than on a uniform grid.

v1: [`packages/treetime-grid/src/interp_nonuniform.rs#L25-L56`](../../packages/treetime-grid/src/interp_nonuniform.rs#L25-L56).

---

## Exponential Damping for Outer-Loop Convergence

Blends optimized branch lengths with previous values using iteration-dependent weights: `bl = bl_new * (1 - damping^(i+1)) + bl_old * damping^(i+1)`. Early iterations take conservative steps; later iterations approach the full Newton update. Prevents oscillation in the alternating optimization (marginal reconstruction / branch length update) cycle.

v1: [`packages/treetime/src/commands/optimize/run.rs`](../../packages/treetime/src/commands/optimize/run.rs) `apply_damping()`.

v0: `optimize_tree_marginal()` at [`packages/legacy/treetime/treetime/treeanc.py#L1297-L1360`](../../packages/legacy/treetime/treetime/treeanc.py#L1297-L1360) with `damping=0.75` default.

References:

- Sagulenko, Puller, and Neher (2018). "TreeTime: Maximum-Likelihood Phylodynamic Analysis." Virus Evolution 4(1):vex042. doi:10.1093/ve/vex042

---

## Zero-Length Branch Detection (Derivative Sign)

Determines if zero is the optimal branch length by evaluating the sign of `d/dt log L(0)`. For independent sites with eigendecomposition-based likelihood `L_i(t) = sum_c k_{ic} exp(lambda_c t)`, the per-site derivative at zero is `(sum_c k_{ic} lambda_c) / (sum_c k_{ic})`. If the total derivative (summed over sites) is negative, the likelihood decreases as `t` increases from zero, making zero a local maximum.

v1: [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs) `is_zero_branch_optimal()`.

v0 does not have an equivalent analytical check. v0 uses `prune_short_branches()` with a compound threshold-and-probability criterion instead.

References:

- Zhang, Dinh, and Matsen (2018). "Non-bifurcating Phylogenetic Tree Inference via the Adaptive LASSO." arXiv:1805.11073 (theoretical background on zero-length branch identification)

---

## Zero-Length Branch Pruning

Collapses internal edges whose optimal length is zero or near-zero, reparenting children to the grandparent. Creates polytomies. v0 uses a compound criterion: `branch_length < 0.1 * one_mutation AND prob_t(parent_seq, child_seq, 0) > 0.1`.

v0: `prune_short_branches()` at [`packages/legacy/treetime/treetime/treeanc.py#L1475-L1496`](../../packages/legacy/treetime/treetime/treeanc.py#L1475-L1496).

v1: Not implemented in the optimize loop. The prune command has `--prune-short` and `--prune-empty` as standalone operations. See known issue `M-optimize-no-topology-cleanup-in-loop.md`.

---

## Shared-Mutation Polytomy Merging

Scans children of polytomy nodes for shared substitutions. When two siblings carry identical mutations, they are grouped under a new internal node with `branch_length = #shared_mutations / alignment_length`. Shared mutations move to the new parent edge; remaining unique mutations stay on child edges.

v1: [`packages/treetime/src/commands/prune/run.rs`](../../packages/treetime/src/commands/prune/run.rs) `merge_shared_mutation_branches()`. Available as `--merge-shared-mutations` on the prune command. Not integrated into the optimize loop.

v0: No formal implementation. Design doc describes "ad-hoc scripts" in nextstrain pathogen pipelines.

---

## Greedy Temporal Polytomy Resolution

For each polytomy, computes pairwise likelihood gain from merging children under a new intermediate node. Uses Brent optimization over the time domain with `zero_branch_slope = mu * L` penalty for the new zero-mutation branch. Greedily picks the best pair above `resolution_threshold` (0.05). O(n^2) per polytomy.

v1: [`packages/treetime/src/commands/timetree/optimization/polytomy.rs`](../../packages/treetime/src/commands/timetree/optimization/polytomy.rs) `resolve_polytomies()`.

v0: `_poly()` at [`packages/legacy/treetime/treetime/treetime.py#L713-L870`](../../packages/legacy/treetime/treetime/treetime.py#L713-L870). v0 distinguishes "stretched" (`mutation_length < clock_length`) from "compressed" children and by default only resolves stretched ones.

---

## Stochastic Coalescent Polytomy Resolution

Simulates a backward-in-time coalescent process. Branches without mutations are "ready to coalesce"; branches with mutations must have mutations removed stochastically. Randomly pairs ready branches for merging. Recommended for large polytomies where greedy mode produces caterpillar-like subtrees.

v0: `generate_subtree()` at [`packages/legacy/treetime/treetime/treetime.py#L872-L1010`](../../packages/legacy/treetime/treetime/treetime.py#L872-L1010).

v1: Not implemented. Tracked: `N-timetree-stochastic-polytomy-unimplemented.md`.

---

## File Index

| File                                                                                                                   | Algorithms                                                         |
| ---------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------ |
| [`packages/treetime/src/commands/optimize/`](../../packages/treetime/src/commands/optimize/)                           | Newton-Raphson, grid search, likelihood eval, damping, zero-detect |
| [`packages/treetime/src/commands/prune/`](../../packages/treetime/src/commands/prune/)                                 | Shared-mutation merging, edge collapsing                           |
| [`packages/treetime/src/commands/timetree/optimization/`](../../packages/treetime/src/commands/timetree/optimization/) | Greedy temporal polytomy resolution                                |
| [`packages/treetime-grid/src/`](../../packages/treetime-grid/src/)                                                     | Interpolation (uniform, non-uniform)                               |
