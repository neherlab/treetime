# Chapter 7: Audit - inventory, inconsistencies, and proposals

[Back to index](_index.md) | Previous: [Chapter 6: The argmin crate](6-argmin-crate.md) | Next: [Glossary](glossary.md)

This chapter is the actionable output of the optimization methods audit. It inventories every optimization in v1, identifies inconsistencies across commands, and proposes prioritized refactoring.

The full item-by-item audit with code locations, callers, callees, and convergence criteria was produced during the initial audit session (no longer available).

## Inventory

### Existing argmin usage (E1-E6)

| ID  | Location                                                                                                                          | Solver                    | Objective              |
| :-- | :-------------------------------------------------------------------------------------------------------------------------------- | :------------------------ | :--------------------- |
| E1  | [method_brent.rs#L36-L78](../../../packages/treetime/src/commands/clock/find_best_root/method_brent.rs#L36-L78)                   | `BrentOpt`                | Root-split chi-squared |
| E2  | [method_golden_section.rs#L36-L81](../../../packages/treetime/src/commands/clock/find_best_root/method_golden_section.rs#L36-L81) | `GoldenSectionSearch`     | Root-split chi-squared |
| E3  | [method_grid_search.rs#L10-L58](../../../packages/treetime/src/commands/clock/find_best_root/method_grid_search.rs#L10-L58)       | Grid search (hand-rolled) | Root-split chi-squared |
| E4  | [optimize_tc.rs#L47-L87](../../../packages/treetime/src/commands/timetree/coalescent/optimize_tc.rs#L47-L87)                      | `BrentOpt`                | Coalescent Tc          |
| E5  | [skyline.rs#L68-L157](../../../packages/treetime/src/commands/timetree/coalescent/skyline.rs#L68-L157)                            | `NelderMead`              | Skyline Tc(t)          |
| E6  | [polytomy.rs#L239-L289](../../../packages/treetime/src/commands/timetree/optimization/polytomy.rs#L239-L289)                      | `BrentOpt`                | Polytomy merge time    |

### Hand-rolled optimization (H1-H11)

| ID  | Location                                                                                                                              | Algorithm             | argmin candidate?                |
| :-- | :------------------------------------------------------------------------------------------------------------------------------------ | :-------------------- | :------------------------------- |
| H1  | [optimize_unified.rs#L533-L711](../../../packages/treetime/src/commands/optimize/optimize_unified.rs#L533-L711)                       | Newton-Raphson + grid | **High** - wrap as `Solver` (P1) |
| H2  | [common.rs#L93-L163](../../../packages/treetime/src/gtr/infer_gtr/common.rs#L93-L163)                                                 | Fixed-point iteration | No - domain-specific ECM         |
| H3  | [run.rs#L275-L410](../../../packages/treetime/src/commands/optimize/run.rs#L275-L410)                                                 | EM-like outer loop    | No - pipeline coordination       |
| H4  | [run.rs#L236-L279](../../../packages/treetime/src/commands/timetree/run.rs#L236-L279)                                                 | Refinement pipeline   | No - multi-phase pipeline        |
| H5  | [distribution.rs#L425-L491](../../../packages/treetime-distribution/src/distribution_core/distribution.rs#L425-L491)                  | HPD bisection         | **Medium** - `BrentRoot` (P2)    |
| H6  | Same as E3                                                                                                                            | Grid search           | Retire (P4)                      |
| H7  | [optimize_unified.rs#L399-L421](../../../packages/treetime/src/commands/optimize/optimize_unified.rs#L399-L421)                       | Grid fallback         | **Medium** - `BrentOpt` (P5)     |
| H8  | [function.rs#L250](../../../packages/treetime-distribution/src/distribution_core/function.rs#L250)                                    | Array argmax          | No - array scan                  |
| H9  | [distribution.rs#L236-L293](../../../packages/treetime-distribution/src/distribution_core/distribution.rs#L236-L293)                  | CDF linear scan       | No - sorted lookup               |
| H10 | [branch_length_likelihood.rs#L31-L63](../../../packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L31-L63) | Grid evaluation       | No - full distribution needed    |
| H11 | [relaxed_clock.rs#L25-L125](../../../packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs#L25-L125)                  | Analytical quadratic  | No - closed-form solution        |

## Inconsistencies

### I1. Three solvers for one problem (root-split)

E1/E2/E3 solve the same problem with different methods. E3 bypasses argmin and calls `evaluate_clock_set()` directly. **Retire E3** in favor of E1/E2.

### I2. Branch length optimizer not on argmin

The most frequently called optimization (H1) is the only major 1D optimizer not using argmin. All other 1D optimizations use `BrentOpt`. **Wrap H1 as argmin `Solver`** (P1).

### I3. Three duplicate observer implementations

E1, E2, E4 each define a near-identical `Observer` struct that logs every 10th iteration. **Extract a shared `PeriodicDebugObserver`** (P3).

### I4. Hardcoded vs configurable parameters

Root-split optimization has full CLI configurability. Coalescent Tc (max_iters=100), polytomy (max_iters=50), and branch length (max_iter=10) are hardcoded. **Make convergence params configurable** (P7).

### I5. Non-concave fallback divergence

H1 uses grid search (100 points) as fallback. All other optimizations use bracket-based methods that handle non-concavity natively. **Replace grid fallback with `BrentOpt`** (P5).

### I6. Skyline: NelderMead (v1) vs SLSQP (v0)

NelderMead has no convergence guarantee above 2D. The skyline cost function has analytically tractable gradients. **Upgrade to `LBFGS`** (P8).

### I7. `&CostFunction` vs `CostFunction` trait pattern

E1/E2/E4/E5 use `impl CostFunction for &Type`. E6 uses `impl CostFunction for Type`. **Standardize** on one pattern.

## Proposals (prioritized)

### P3. Extract generic observer (zero risk)

Eliminate 3 near-identical observer structs. Pure deduplication.

### P5. BrentOpt fallback (low risk)

Replace the 100-point grid search fallback in Newton branch optimization with `BrentOpt` on `[0.1 * one_mutation, 1.5 * branch_length + one_mutation]`. Guaranteed convergence, higher precision, proven pattern in this codebase.

Risk: if the surface has multiple local minima (possible for GTR models per <a id="cite-1"></a>[Dinh and Matsen 2017](https://doi.org/10.1214/16-AAP1240) [[1](#ref-1)]), BrentOpt could find a different local minimum than the grid. Mitigation: use the grid result as initial bracket center.

### P4. Retire grid search (zero risk)

Remove `method_grid_search.rs` and the `Grid` variant from `BranchPointOptimizationParams`. BrentOpt and GoldenSection solve the same problem more efficiently.

### P1. Newton as argmin Solver (low risk)

Wrap the branch length Newton-Raphson as a custom `argmin::core::Solver`. The mathematical computation stays identical - only the iteration scaffolding changes. Enables swapping to BrentOpt without changing callers, and gets argmin's observer and convergence infrastructure.

### P2. BrentRoot for HPD (low risk)

Replace hand-rolled bisection in HPD computation with argmin `BrentRoot`. Primarily for consistency - the current bisection works and is not a performance bottleneck.

### P8. Skyline LBFGS (medium risk)

Replace NelderMead with LBFGS for skyline optimization. Requires deriving and validating analytical gradients of the coalescent likelihood + GMRF penalty + boundary penalty. The gradient is tractable (see [Chapter 3](3-coalescent-skyline.md)) but incorrect gradients cause LBFGS to diverge.

### P6. Shared optimization module (depends on P1/P3)

Create a shared module for generic observer, cost function patterns, optimization parameter types, and argmin utilities.

### P7. Configurable convergence parameters

Add CLI flags for Tc, branch-length, and polytomy convergence parameters with current hardcoded values as defaults.

## TreeTime v0 optimization inventory

Complete `scipy.optimize` call sites in v0 for comparison:

| v0 location             | scipy method                | Purpose             | v1 equivalent       |
| :---------------------- | :-------------------------- | :------------------ | :------------------ |
| `gtr.py:879`            | `minimize_scalar` (brent)   | Branch length       | Newton-Raphson (H1) |
| `merger_models.py:260`  | `minimize_scalar` (brent)   | Constant Tc         | BrentOpt (E4)       |
| `merger_models.py:310`  | `minimize` (SLSQP)          | Skyline Tc(t)       | NelderMead (E5)     |
| `treeregression.py:340` | `minimize_scalar` (bounded) | Root position       | BrentOpt (E1)       |
| `clock_tree.py:1187`    | `minimize_scalar` (brent)   | HPD interval        | Bisection (H5)      |
| `treetime.py:748`       | `minimize_scalar` (bounded) | Polytomy merge time | BrentOpt (E6)       |
| `treeanc.py:1683`       | `minimize_scalar` (brent)   | GTR rate mu         | Not ported          |

v0 features not in v1 (optimization-related): `optimize_gtr_rate()`, `infer_gtr_iterative()`, `generate_subtree()` (stochastic polytomy), incremental gain matrix update in polytomy resolution, stretched vs compressed branch separation, skyline confidence via numerical second derivatives.

## References

1. <a id="ref-1"></a> Dinh, Vu, and Frederick A. Matsen IV. 2017. "The Shape of the One-Dimensional Phylogenetic Likelihood Function." _Ann. Appl. Prob._ 27(3):1646-1677. https://doi.org/10.1214/16-AAP1240 [↩](#cite-1)
