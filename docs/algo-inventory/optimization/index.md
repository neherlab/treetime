# Numerical Optimization Algorithms

[Back to index](../)

## Newton-Raphson

| Property    | Value                                                                  |
| ----------- | ---------------------------------------------------------------------- |
| Type        | Well-known                                                             |
| v1 Location | `packages/treetime/src/commands/optimize/optimize_unified.rs:231-251:` |
| Reference   | Nocedal & Wright, "Numerical Optimization." Chapter 2                  |
| Reference   | Felsenstein (2003). "Inferring Phylogenies." Chapter 16                |

`t_new = t - clamp(f'/f'', -1.0, t)` with fallback to grid search.

---

## Eigendecomposition-Based Likelihood

| Property    | Value                                                                                       |
| ----------- | ------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                  |
| v1 Location | `packages/treetime/src/commands/optimize/optimize_dense_eval.rs`, `optimize_sparse_eval.rs` |
| Reference   | Felsenstein (1981). "Evolutionary trees from DNA sequences."                                |
| Reference   | Yang (2006). "Computational Molecular Evolution." Chapter 4                                 |

Precomputes eigenvector coefficients for efficient likelihood + derivative evaluation.

---

## Piecewise Linear Interpolation (Uniform)

| Property    | Value                                            |
| ----------- | ------------------------------------------------ |
| Type        | Well-known                                       |
| v1 Location | `packages/treetime-grid/src/grid_fn.rs:279-351:` |
| Reference   | Burden & Faires, "Numerical Analysis." Chapter 3 |

O(1) interval lookup via `floor((x - x_min) / dx)`.

---

## Piecewise Linear Interpolation (Non-Uniform)

| Property    | Value                                                    |
| ----------- | -------------------------------------------------------- |
| Type        | Well-known                                               |
| v1 Location | `packages/treetime-grid/src/interp_nonuniform.rs:25-56:` |

O(log n) binary search interval lookup.

---

## File Index

| File                                           | Algorithms                                   |
| ---------------------------------------------- | -------------------------------------------- |
| `packages/treetime/src/commands/optimize/*.rs` | Newton-Raphson, grid search, likelihood eval |
| `packages/treetime-grid/src/*.rs`              | Interpolation (uniform, non-uniform)         |
