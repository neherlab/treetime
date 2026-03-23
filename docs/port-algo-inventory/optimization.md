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

## Piecewise Linear Interpolation (Uniform Grid)

O(1) interval lookup via `floor((x - x_min) / dx)` for uniformly spaced grids. Used throughout the distribution system for evaluating discretized probability distributions and branch length likelihoods on fixed grids.

v1: [`packages/treetime-grid/src/grid_fn.rs#L279-L351`](../../packages/treetime-grid/src/grid_fn.rs#L279-L351).

Reference: Burden & Faires. "Numerical Analysis." Chapter 3.

---

## Piecewise Linear Interpolation (Non-Uniform Grid)

O(log n) binary search interval lookup for non-uniformly spaced grids. Used for the skyline coalescent Tc(t) function where grid points are placed at coalescent event times rather than on a uniform grid.

v1: [`packages/treetime-grid/src/interp_nonuniform.rs#L25-L56`](../../packages/treetime-grid/src/interp_nonuniform.rs#L25-L56).

---

## File Index

| File                                                                                         | Algorithms                                   |
| -------------------------------------------------------------------------------------------- | -------------------------------------------- |
| [`packages/treetime/src/commands/optimize/`](../../packages/treetime/src/commands/optimize/) | Newton-Raphson, grid search, likelihood eval |
| [`packages/treetime-grid/src/`](../../packages/treetime-grid/src/)                           | Interpolation (uniform, non-uniform)         |
