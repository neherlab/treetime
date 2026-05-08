# Clock Inference Algorithms

[Back to index](README.md)

## WLS Sufficient Statistics

TreeTime estimates the molecular clock rate via weighted least squares (WLS) regression of root-to-tip divergence against sampling dates (<a id="cite-1"></a>[Sagulenko, Puller, and Neher 2018](https://doi.org/10.1093/ve/vex042) [[1](#ref-1)], Equations 11-14). Rather than materializing the full N x N covariance matrix for generalized least squares, TreeTime propagates six sufficient statistics through the tree in O(N) time using the `ClockSet` data structure:

- `t_sum`: sum of weighted dates
- `tsq_sum`: sum of weighted squared dates
- `d_sum`: sum of weighted divergences
- `dsq_sum`: sum of weighted squared divergences
- `dt_sum`: sum of weighted date-divergence cross products
- `norm`: sum of weights (effective sample size)

From these six values at any node, the clock rate and intercept follow directly: `rate = (norm * dt_sum - t_sum * d_sum) / (norm * tsq_sum - t_sum^2)`. The chi-squared statistic, R-value, and parameter covariance matrix (2x2 Hessian inverse) are also derivable from the sufficient statistics without additional tree traversal.

v1: [`packages/treetime/src/commands/clock/clock_set.rs#L1-L172`](../../packages/treetime/src/commands/clock/clock_set.rs#L1-L172).
v0: [`packages/legacy/treetime/treetime/treeregression.py#L7-L300`](../../packages/legacy/treetime/treetime/treeregression.py#L7-L300).

Key functions: `ClockSet::clock_rate()` (`#ClockSet`, `#clock_rate`), `ClockSet::chisq()` (`#chisq`), `ClockSet::propagate_averages()` (`#propagate_averages`).

### References

- <a id="ref-1"></a>Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4(1):vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-1)

---

## Tree Regression

Two-pass message passing that implements phylogenetic GLS regression (<a id="cite-2"></a>[Felsenstein 1985](https://doi.org/10.1086/284325) [[2](#ref-2)]) in O(N) time. Standard OLS assumes independent residuals, but in phylogenetic data, samples sharing internal branches have correlated root-to-tip distances. The tree regression accounts for this by propagating weighted statistics along the tree structure, equivalent to GLS with a phylogenetic covariance matrix but without materializing the N x N matrix.

The backward pass (postorder) accumulates `ClockSet` statistics from leaves toward the root. Each leaf contributes its date and divergence. Internal nodes merge children's statistics via `propagate_averages()`. The forward pass (preorder) propagates root-level statistics back down, enabling the best root position search: at each node, the regression can be re-evaluated as if that node were the root by combining the subtree statistics (from the backward pass) with the complement statistics (from the forward pass).

v1: [`packages/treetime/src/commands/clock/clock_regression.rs#L47-L125`](../../packages/treetime/src/commands/clock/clock_regression.rs#L47-L125).
v0: [`packages/legacy/treetime/treetime/treeregression.py#L192-L283`](../../packages/legacy/treetime/treetime/treeregression.py#L192-L283).

Re-estimation mode: when `prev_clock_rate` is provided (during iterative timetree refinement), the regression uses solver-updated time lengths converted to divergence (`time_length * rate * gamma`) instead of input branch lengths. The `edge_divergence()` helper (`#edge_divergence`, [`packages/treetime/src/commands/clock/clock_regression.rs#L238-L250`](../../packages/treetime/src/commands/clock/clock_regression.rs#L238-L250)) handles this fallback logic.

### Variance model

The `ClockParams` (`#ClockParams`) struct at [`packages/treetime/src/commands/clock/clock_regression.rs#L18-L34`](../../packages/treetime/src/commands/clock/clock_regression.rs#L18-L34) controls the variance model: `variance_factor` (branch-length-proportional term), `variance_offset` (constant term for internal branches), `variance_offset_leaf` (additional offset for terminal branches, controlled by `--tip-slack`). When `--covariation` is active, these parameters are set to account for phylogenetic correlation in the regression.

Key functions: `clock_regression_backward()` (`#clock_regression_backward`), `clock_regression_forward()` (`#clock_regression_forward`), `estimate_clock_model_with_reroot()` (`#estimate_clock_model_with_reroot`), `estimate_clock_model_with_reroot_policy()` (`#estimate_clock_model_with_reroot_policy`).

### References

- <a id="ref-2"></a>Felsenstein, Joseph. 1985. "Phylogenies and the Comparative Method." _American Naturalist_ 125(1):1-15. https://doi.org/10.1086/284325 [↩](#cite-2)
- Sagulenko, Puller, and Neher 2018 [[1](#ref-1)]

---

## Brent's Method

Hybrid 1D optimization (<a id="cite-3"></a>[Brent 1973](https://doi.org/10.1007/978-3-0348-5952-3) [[3](#ref-3)]) combining parabolic interpolation (fast convergence near minimum) with golden section search (guaranteed convergence). Used for optimizing the root position along a branch: given a branch with endpoints, Brent's method finds the split point that minimizes the clock regression chi-squared.

v1: [`packages/treetime/src/commands/clock/find_best_root/method_brent.rs#L36-L78`](../../packages/treetime/src/commands/clock/find_best_root/method_brent.rs#L36-L78).
v0: `scipy.optimize.minimize_scalar` with `method='bounded'`.

### References

- <a id="ref-3"></a>Brent, Richard P. 1973. _Algorithms for Minimization Without Derivatives._ Prentice-Hall. ISBN 978-0-13-022335-7. [↩](#cite-3)

---

## Golden Section Search

Bracket-based 1D optimization (<a id="cite-4"></a>[Kiefer 1953](https://doi.org/10.2307/2032012) [[4](#ref-4)]) with O(log(1/epsilon)) convergence. At each step, the bracket is narrowed by the golden ratio phi = (1+sqrt(5))/2, maintaining the ratio between the two sub-intervals. Slower than Brent's method (which uses parabolic acceleration) but simpler and guaranteed to converge without requiring derivative information.

v1: [`packages/treetime/src/commands/clock/find_best_root/method_golden_section.rs#L36-L81`](../../packages/treetime/src/commands/clock/find_best_root/method_golden_section.rs#L36-L81).

### References

- <a id="ref-4"></a>Kiefer, Jack. 1953. "Sequential Minimax Search for a Maximum." _Proceedings of the American Mathematical Society_ 4(3):502-506. https://doi.org/10.2307/2032012 [↩](#cite-4)

---

## IQD Outlier Detection

Identifies clock outliers using Tukey's fences (<a id="cite-5"></a>[Tukey 1977](https://doi.org/10.1007/978-1-4419-7976-6) [[5](#ref-5)]) applied to clock deviation residuals. The interquartile distance IQD = Q3 - Q1 measures the spread of the middle 50% of clock deviations. A leaf with `|clock_deviation| > threshold * IQD` is flagged as an outlier, where the threshold is controlled by `--clock-filter` (default 3.0).

IQD-based detection is nonparametric and does not assume normally distributed residuals, making it appropriate for phylogenetic data where residual distributions can be skewed. For comparison, under a normal distribution IQD = 1.35 _ sigma, so `threshold _ IQD = 3 _ 1.35 _ sigma = 4.05 \* sigma` - roughly the 99.995% interval.

v1: [`packages/treetime/src/commands/clock/clock_filter.rs#L21-L105`](../../packages/treetime/src/commands/clock/clock_filter.rs#L21-L105).
v0: [`packages/legacy/treetime/treetime/clock_filter_methods.py#L5-L40`](../../packages/legacy/treetime/treetime/clock_filter_methods.py#L5-L40).

### References

- <a id="ref-5"></a>Tukey, John W. 1977. _Exploratory Data Analysis._ Addison-Wesley. ISBN 978-0-201-07616-5. [↩](#cite-5)

---

## Tree Rerooting

Edge inversion along the path from the old root to the new root position, with cleanup of trivial (single-child) nodes. Rerooting changes which node is the root of the tree without altering the unrooted topology. All edges along the old-root-to-new-root path are inverted (parent becomes child and vice versa), and the old root is removed if it becomes a trivial single-child node.

v1: [`packages/treetime/src/commands/clock/reroot.rs#L82-L283`](../../packages/treetime/src/commands/clock/reroot.rs#L82-L283).

Key functions: `reroot_in_place()` (`#reroot_in_place`), `apply_reroot()` (`#apply_reroot`), `create_new_root_node()` (`#create_new_root_node`).

---

## Best Root Search

Scans all nodes to find the root position that maximizes temporal signal (positive clock rate, minimal chi-squared residuals from the WLS regression). For each candidate node, the regression is evaluated using the pre-computed forward-pass statistics (no re-traversal needed). Only candidates with positive clock rates are considered; negative rates indicate the tree cannot be clock-like when rooted at that position.

After finding the best node, the algorithm refines the root position along the parent and child branches via `find_best_split()`, which uses Brent's method or golden section search to optimize the exact split point along the edge.

v1: [`packages/treetime/src/commands/clock/find_best_root/find_best_root.rs#L19-L135`](../../packages/treetime/src/commands/clock/find_best_root/find_best_root.rs#L19-L135).

Key functions: `find_best_root()` (`#find_best_root`), `has_positive_clock_rate()` (`#has_positive_clock_rate`).

---

## Additional Algorithms

- Grid Search: [`packages/treetime/src/commands/clock/find_best_root/method_grid_search.rs`](../../packages/treetime/src/commands/clock/find_best_root/method_grid_search.rs) - Brute-force 1D optimization over a uniform grid of split positions. Fallback when bracket-based methods cannot be applied.
- Branch Point Cost Function: [`packages/treetime/src/commands/clock/find_best_root/cost_function.rs`](../../packages/treetime/src/commands/clock/find_best_root/cost_function.rs) - Chi-squared objective function for root position optimization, computing regression residuals at a candidate split point.
- Root-to-Tip Distance: [`packages/treetime/src/commands/clock/rtt.rs`](../../packages/treetime/src/commands/clock/rtt.rs) - Tree traversal accumulating branch lengths from root to each leaf, producing the divergence values used in clock regression.
- Clock Model: [`packages/treetime/src/commands/clock/clock_model.rs`](../../packages/treetime/src/commands/clock/clock_model.rs) - Rate/intercept estimation with `ClockModelStats::Estimated` (from regression) or `ClockModelStats::Fixed` (from `--clock-rate`) variants.
- Clock Traits: [`packages/treetime/src/commands/clock/clock_traits.rs`](../../packages/treetime/src/commands/clock/clock_traits.rs) - `ClockNode` and `ClockEdge` trait definitions. `ClockEdge::gamma()` provides per-branch relaxed clock rate multiplier (default 1.0).
- Date Constraints: [`packages/treetime/src/commands/clock/date_constraints.rs`](../../packages/treetime/src/commands/clock/date_constraints.rs) - `load_date_constraints()` (`#load_date_constraints`) assigns date distributions to tree nodes with validation.

---

## Unimplemented

See [unimplemented](unimplemented.md) for full details:

- Local outlier filter (z-score based)
- Full covariance matrix computation (TreeRegression)
- Numerical Hessian for root position uncertainty

---

## File Index

| File                                                                                                                                                   | Algorithms                                               |
| ------------------------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------- |
| [`packages/treetime/src/commands/clock/clock_set.rs`](../../packages/treetime/src/commands/clock/clock_set.rs)                                         | Sufficient statistics WLS                                |
| [`packages/treetime/src/commands/clock/clock_regression.rs`](../../packages/treetime/src/commands/clock/clock_regression.rs)                           | Tree message passing, edge divergence, reroot estimation |
| [`packages/treetime/src/commands/clock/clock_model.rs`](../../packages/treetime/src/commands/clock/clock_model.rs)                                     | Clock model estimation                                   |
| [`packages/treetime/src/commands/clock/clock_traits.rs`](../../packages/treetime/src/commands/clock/clock_traits.rs)                                   | ClockNode, ClockEdge traits                              |
| [`packages/treetime/src/commands/clock/find_best_root/`](../../packages/treetime/src/commands/clock/find_best_root/)                                   | Best root search, Brent, golden section, grid search     |
| [`packages/treetime/src/commands/clock/find_best_root/find_best_root.rs`](../../packages/treetime/src/commands/clock/find_best_root/find_best_root.rs) | Best root node selection                                 |
| [`packages/treetime/src/commands/clock/clock_filter.rs`](../../packages/treetime/src/commands/clock/clock_filter.rs)                                   | IQD outlier detection                                    |
| [`packages/treetime/src/commands/clock/reroot.rs`](../../packages/treetime/src/commands/clock/reroot.rs)                                               | Tree rerooting                                           |
| [`packages/treetime/src/commands/clock/date_constraints.rs`](../../packages/treetime/src/commands/clock/date_constraints.rs)                           | Date constraint loading and validation                   |
| [`packages/treetime/src/commands/clock/rtt.rs`](../../packages/treetime/src/commands/clock/rtt.rs)                                                     | Root-to-tip distance computation                         |
| [`packages/treetime/src/commands/clock/run.rs`](../../packages/treetime/src/commands/clock/run.rs)                                                     | Clock command orchestration                              |
