# Clock Inference Algorithms

[Back to index](../index.md)

## WLS Sufficient Stats

| Property    | Value                                                                                                                                           |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Custom (TreeTime-specific)                                                                                                                      |
| v1 Location | [`packages/treetime/src/commands/clock/clock_set.rs#L1-L172`](../../../packages/treetime/src/commands/clock/clock_set.rs#L1-L172)               |
| v0 Location | [`packages/legacy/treetime/treetime/treeregression.py#L7-L300`](../../../packages/legacy/treetime/treetime/treeregression.py#L7-L300)           |
| Functions   | `ClockSet::clock_rate()` (`#ClockSet`, `#clock_rate`), `ClockSet::chisq()` (`#chisq`), `ClockSet::propagate_averages()` (`#propagate_averages`) |
| Reference   | Neher et al. (2018). "TreeTime: Maximum-likelihood phylodynamic analysis." Virus Evolution, 4(1):vex042, Equations 11-14                        |
| Paper URL   | https://doi.org/10.1093/ve/vex042                                                                                                               |

Propagates 6 sufficient statistics (t_sum, tsq_sum, d_sum, dsq_sum, dt_sum, norm) through tree for O(N) phylogenetic regression.

---

## Tree Regression

| Property    | Value                                                                                                                                                                                                                                                                                          |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known (belief propagation)                                                                                                                                                                                                                                                                |
| v1 Location | [`packages/treetime/src/commands/clock/clock_regression.rs#L47-L125`](../../../packages/treetime/src/commands/clock/clock_regression.rs#L47-L125)                                                                                                                                              |
| v0 Location | [`packages/legacy/treetime/treetime/treeregression.py#L192-L283`](../../../packages/legacy/treetime/treetime/treeregression.py#L192-L283)                                                                                                                                                      |
| Functions   | `clock_regression_backward()` (`#clock_regression_backward`), `clock_regression_forward()` (`#clock_regression_forward`), `estimate_clock_model_with_reroot()` (`#estimate_clock_model_with_reroot`), `estimate_clock_model_with_reroot_policy()` (`#estimate_clock_model_with_reroot_policy`) |
| Reference   | Felsenstein (1985). "Phylogenies and the comparative method." American Naturalist, 125(1):1-15; Neher 2018                                                                                                                                                                                     |

Two-pass message passing implementing GLS regression in O(N) time. Supports re-estimation mode via `prev_clock_rate` parameter: when provided, regression uses solver-updated time lengths converted to divergence (`time_length * rate * gamma`) instead of input branch lengths. The `edge_divergence()` helper (`#edge_divergence`, [L238-L250](../../../packages/treetime/src/commands/clock/clock_regression.rs#L238-L250)) handles this fallback logic.

---

## Brent's Method

| Property    | Value                                                                                                                                                                 |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                                                            |
| v1 Location | [`packages/treetime/src/commands/clock/find_best_root/method_brent.rs#L36-L78`](../../../packages/treetime/src/commands/clock/find_best_root/method_brent.rs#L36-L78) |
| v0 Location | scipy `minimize_scalar` with `method='bounded'`                                                                                                                       |
| Reference   | Brent, R.P. (1973). "Algorithms for Minimization Without Derivatives." Prentice-Hall, Chapter 5                                                                       |

Hybrid parabolic interpolation + golden section search.

---

## Golden Section Search

| Property    | Value                                                                                                                                                                                   |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                                                                              |
| v1 Location | [`packages/treetime/src/commands/clock/find_best_root/method_golden_section.rs#L36-L81`](../../../packages/treetime/src/commands/clock/find_best_root/method_golden_section.rs#L36-L81) |
| Reference   | Kiefer, J. (1953). "Sequential minimax search for a maximum." Proc AMS, 4(3):502-506                                                                                                    |

Bracket-based 1D optimization with O(log(1/epsilon)) convergence.

---

## IQD Outlier Detection

| Property    | Value                                                                                                                                           |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known (Tukey's fences)                                                                                                                     |
| v1 Location | [`packages/treetime/src/commands/clock/clock_filter.rs#L21-L105`](../../../packages/treetime/src/commands/clock/clock_filter.rs#L21-L105)       |
| v0 Location | [`packages/legacy/treetime/treetime/clock_filter_methods.py#L5-L40`](../../../packages/legacy/treetime/treetime/clock_filter_methods.py#L5-L40) |
| Reference   | Tukey, J.W. (1977). "Exploratory Data Analysis." Addison-Wesley                                                                                 |

Flags outliers where `|clock_deviation| > threshold * IQD`.

---

## Tree Rerooting

| Property    | Value                                                                                                                            |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Standard graph operation                                                                                                         |
| v1 Location | [`packages/treetime/src/commands/clock/reroot.rs#L82-L283`](../../../packages/treetime/src/commands/clock/reroot.rs#L82-L283)    |
| Functions   | `reroot_in_place()` (`#reroot_in_place`), `apply_reroot()` (`#apply_reroot`), `create_new_root_node()` (`#create_new_root_node`) |

Edge inversion along path from old root to new root, with trivial node cleanup.

---

## Best Root Search

| Property    | Value                                                                                                                                                                       |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Custom (TreeTime-specific)                                                                                                                                                  |
| v1 Location | [`packages/treetime/src/commands/clock/find_best_root/find_best_root.rs#L19-L135`](../../../packages/treetime/src/commands/clock/find_best_root/find_best_root.rs#L19-L135) |
| Functions   | `find_best_root()` (`#find_best_root`), `has_positive_clock_rate()` (`#has_positive_clock_rate`)                                                                            |

Iterates all nodes to find the one with lowest chi-squared and positive clock rate, then optimizes position along parent and child branches via `find_best_split()`. Returns an actionable error when all root positions yield a negative clock rate.

---

## Additional Algorithms

- **Grid Search** ([`packages/treetime/src/commands/clock/find_best_root/method_grid_search.rs`](../../../packages/treetime/src/commands/clock/find_best_root/method_grid_search.rs)): Brute-force 1D optimization
- **Branch Point Cost Function** ([`packages/treetime/src/commands/clock/find_best_root/cost_function.rs`](../../../packages/treetime/src/commands/clock/find_best_root/cost_function.rs)): Chi-squared objective for root optimization
- **Variance Model** ([`packages/treetime/src/commands/clock/clock_regression.rs#L18-L34`](../../../packages/treetime/src/commands/clock/clock_regression.rs#L18-L34)): Covariation-aware weighting (`ClockParams` (`#ClockParams`): `variance_factor`, `variance_offset`, `variance_offset_leaf`)
- **Root-to-Tip Distance** ([`packages/treetime/src/commands/clock/rtt.rs`](../../../packages/treetime/src/commands/clock/rtt.rs)): Tree traversal accumulation
- **Clock Model** ([`packages/treetime/src/commands/clock/clock_model.rs`](../../../packages/treetime/src/commands/clock/clock_model.rs)): Rate/intercept estimation with `ClockModelStats::Estimated` or `ClockModelStats::Fixed` variants
- **Clock Traits** ([`packages/treetime/src/commands/clock/clock_traits.rs`](../../../packages/treetime/src/commands/clock/clock_traits.rs)): `ClockNode` and `ClockEdge` trait definitions; `ClockEdge::gamma()` provides per-branch relaxed clock rate multiplier (default 1.0)
- **Date Constraints** ([`packages/treetime/src/commands/clock/date_constraints.rs`](../../../packages/treetime/src/commands/clock/date_constraints.rs)): `load_date_constraints()` (`#load_date_constraints`) assigns date distributions to tree nodes with validation

---

## Unimplemented

See [unimplemented](../unimplemented/index.md) for full details:

- Local outlier filter (z-score based)
- Full covariance matrix computation
- Rate susceptibility analysis
- Numerical Hessian for root position uncertainty

---

## File Index

| File                                                                                                                                                      | Algorithms                                               |
| --------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------- |
| [`packages/treetime/src/commands/clock/clock_set.rs`](../../../packages/treetime/src/commands/clock/clock_set.rs)                                         | Sufficient statistics WLS                                |
| [`packages/treetime/src/commands/clock/clock_regression.rs`](../../../packages/treetime/src/commands/clock/clock_regression.rs)                           | Tree message passing, edge divergence, reroot estimation |
| [`packages/treetime/src/commands/clock/clock_model.rs`](../../../packages/treetime/src/commands/clock/clock_model.rs)                                     | Clock model estimation                                   |
| [`packages/treetime/src/commands/clock/clock_traits.rs`](../../../packages/treetime/src/commands/clock/clock_traits.rs)                                   | ClockNode, ClockEdge traits                              |
| [`packages/treetime/src/commands/clock/find_best_root/`](../../../packages/treetime/src/commands/clock/find_best_root/)                                   | Best root search, Brent, golden section, grid search     |
| [`packages/treetime/src/commands/clock/find_best_root/find_best_root.rs`](../../../packages/treetime/src/commands/clock/find_best_root/find_best_root.rs) | Best root node selection                                 |
| [`packages/treetime/src/commands/clock/clock_filter.rs`](../../../packages/treetime/src/commands/clock/clock_filter.rs)                                   | IQD outlier detection                                    |
| [`packages/treetime/src/commands/clock/reroot.rs`](../../../packages/treetime/src/commands/clock/reroot.rs)                                               | Tree rerooting                                           |
| [`packages/treetime/src/commands/clock/date_constraints.rs`](../../../packages/treetime/src/commands/clock/date_constraints.rs)                           | Date constraint loading and validation                   |
| [`packages/treetime/src/commands/clock/rtt.rs`](../../../packages/treetime/src/commands/clock/rtt.rs)                                                     | Root-to-tip distance computation                         |
| [`packages/treetime/src/commands/clock/run.rs`](../../../packages/treetime/src/commands/clock/run.rs)                                                     | Clock command orchestration                              |
