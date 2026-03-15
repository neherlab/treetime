# Numerical Optimization Algorithms

[Back to index](_index.md)

## Newton-Raphson

| Property    | Value                                                                                                                                                  |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Type        | Well-known                                                                                                                                             |
| v1 Location | [`packages/treetime/src/commands/optimize/optimize_unified.rs#L231-L251`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L231-L251) |
| Reference   | Nocedal & Wright, "Numerical Optimization." Chapter 2                                                                                                  |
| Reference   | Felsenstein (2003). "Inferring Phylogenies." Chapter 16                                                                                                |

`t_new = t - clamp(f'/f'', -1.0, t)` with fallback to grid search.

---

## Eigendecomposition-Based Likelihood

| Property    | Value                                                                                                                                                                                                                                                                                |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Type        | Well-known                                                                                                                                                                                                                                                                           |
| v1 Location | [`packages/treetime/src/commands/optimize/optimize_dense_eval.rs`](../../packages/treetime/src/commands/optimize/optimize_dense_eval.rs), [`packages/treetime/src/commands/optimize/optimize_sparse_eval.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse_eval.rs) |
| Reference   | Felsenstein (1981). "Evolutionary trees from DNA sequences."                                                                                                                                                                                                                         |
| Reference   | Yang (2006). "Computational Molecular Evolution." Chapter 4                                                                                                                                                                                                                          |

Precomputes eigenvector coefficients for efficient likelihood + derivative evaluation.

---

## Piecewise Linear Interpolation (Uniform)

| Property    | Value                                                                                                      |
| ----------- | ---------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                 |
| v1 Location | [`packages/treetime-grid/src/grid_fn.rs#L279-L351`](../../packages/treetime-grid/src/grid_fn.rs#L279-L351) |
| Reference   | Burden & Faires, "Numerical Analysis." Chapter 3                                                           |

O(1) interval lookup via `floor((x - x_min) / dx)`.

---

## Piecewise Linear Interpolation (Non-Uniform)

| Property    | Value                                                                                                                      |
| ----------- | -------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                 |
| v1 Location | [`packages/treetime-grid/src/interp_nonuniform.rs#L25-L56`](../../packages/treetime-grid/src/interp_nonuniform.rs#L25-L56) |

O(log n) binary search interval lookup.

---

## File Index

| File                                                                                         | Algorithms                                   |
| -------------------------------------------------------------------------------------------- | -------------------------------------------- |
| [`packages/treetime/src/commands/optimize/`](../../packages/treetime/src/commands/optimize/) | Newton-Raphson, grid search, likelihood eval |
| [`packages/treetime-grid/src/`](../../packages/treetime-grid/src/)                           | Interpolation (uniform, non-uniform)         |
