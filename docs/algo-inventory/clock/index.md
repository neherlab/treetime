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

| Property    | Value                                                                                                                                             |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known (belief propagation)                                                                                                                   |
| v1 Location | [`packages/treetime/src/commands/clock/clock_regression.rs#L42-L113`](../../../packages/treetime/src/commands/clock/clock_regression.rs#L42-L113) |
| v0 Location | [`packages/legacy/treetime/treetime/treeregression.py#L192-L283`](../../../packages/legacy/treetime/treetime/treeregression.py#L192-L283)         |
| Reference   | Felsenstein (1985). "Phylogenies and the comparative method." American Naturalist, 125(1):1-15; Neher 2018                                        |

Two-pass message passing implementing GLS regression in O(N) time.

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
| v1 Location | [`packages/treetime/src/commands/clock/reroot.rs#L82-L282`](../../../packages/treetime/src/commands/clock/reroot.rs#L82-L282)    |
| Functions   | `reroot_in_place()` (`#reroot_in_place`), `apply_reroot()` (`#apply_reroot`), `create_new_root_node()` (`#create_new_root_node`) |

Edge inversion along path from old root to new root, with trivial node cleanup.

---

## Additional Algorithms

- **Grid Search** ([`packages/treetime/src/commands/clock/find_best_root/method_grid_search.rs`](../../../packages/treetime/src/commands/clock/find_best_root/method_grid_search.rs)): Brute-force 1D optimization
- **Branch Point Cost Function** ([`packages/treetime/src/commands/clock/find_best_root/cost_function.rs`](../../../packages/treetime/src/commands/clock/find_best_root/cost_function.rs)): Chi-squared objective for root optimization
- **Variance Model** ([`packages/treetime/src/commands/clock/clock_regression.rs#L17-L33`](../../../packages/treetime/src/commands/clock/clock_regression.rs#L17-L33)): Covariation-aware weighting
- **Root-to-Tip Distance** ([`packages/treetime/src/commands/clock/rtt.rs`](../../../packages/treetime/src/commands/clock/rtt.rs)): Tree traversal accumulation

---

## Unimplemented

See [unimplemented](../unimplemented/index.md) for full details:

- Local outlier filter (z-score based)
- Full covariance matrix computation
- Rate susceptibility analysis
- Numerical Hessian for root position uncertainty

---

## File Index

| File                                                                                                                            | Algorithms                         |
| ------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------- |
| [`packages/treetime/src/commands/clock/clock_set.rs`](../../../packages/treetime/src/commands/clock/clock_set.rs)               | Sufficient statistics WLS          |
| [`packages/treetime/src/commands/clock/clock_regression.rs`](../../../packages/treetime/src/commands/clock/clock_regression.rs) | Tree message passing               |
| [`packages/treetime/src/commands/clock/find_best_root/`](../../../packages/treetime/src/commands/clock/find_best_root/)         | Brent, golden section, grid search |
| [`packages/treetime/src/commands/clock/clock_filter.rs`](../../../packages/treetime/src/commands/clock/clock_filter.rs)         | IQD outlier detection              |
| [`packages/treetime/src/commands/clock/reroot.rs`](../../../packages/treetime/src/commands/clock/reroot.rs)                     | Tree rerooting                     |
