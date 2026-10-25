# Clock Inference Algorithms

## Least-Squares Sufficient Statistics

TreeTime estimates the molecular clock rate by an ordinary least-squares (OLS) regression of root-to-tip divergence against sampling date (<a id="cite-1"></a>[Sagulenko, Puller, and Neher 2018](https://doi.org/10.1093/ve/vex042) [[1](#ref-1)], Equations 11-14), with one unit weight per dated tip that is not an outlier ([kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md)). The `ClockSet` data structure holds six sums over those tips:

- `t_sum`: sum of dates
- `tsq_sum`: sum of squared dates
- `d_sum`: sum of divergences
- `dsq_sum`: sum of squared divergences
- `dt_sum`: sum of date-divergence products
- `norm`: number of tips

From these six values, the clock rate and intercept follow directly: `rate = (norm * dt_sum - t_sum * d_sum) / (norm * tsq_sum - t_sum^2)`; the chi-squared statistic (half the residual sum of squares) and the R-value need no further traversal. Moving the root by $\ell$ lengthens the divergences of all tips in a set by $\ell$, which `ClockSet::shift()` applies to the sums: $\sum d \mathrel{+}= \ell n$, $\sum td \mathrel{+}= \ell \sum t$, $\sum d^2 \mathrel{+}= 2\ell\sum d + \ell^2 n$.

v1: `struct ClockSet` [packages/treetime/src/clock/clock_set.rs](../../packages/treetime/src/clock/clock_set.rs).
v0: [`packages/legacy/treetime/treetime/treeregression.py#L7-L300`](../../packages/legacy/treetime/treetime/treeregression.py#L7-L300).

Key functions: `ClockSet::clock_rate()`, `ClockSet::chisq()`, `ClockSet::chisq_fixed_rate()`, `ClockSet::shift()`.

The regression reports no covariance: OLS treats the tips as independent while they share ancestry. The rate uncertainty comes from the Gaussian timetree at the final rate (`fn gaussian_clock()`, [packages/treetime/src/clock/gaussian.rs](../../packages/treetime/src/clock/gaussian.rs)) and is reported as `rate_std` of an estimated rate when the sequence length is known.

### References

- <a id="ref-1"></a>Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4(1):vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-1)

---

## Tree Regression

Two-pass message passing gives the regression sums for every root position in O(N) time. The backward pass (postorder) accumulates `ClockSet` sums from the leaves: each dated, non-outlier leaf contributes its date, and every node passes its sums to its parent shifted by the branch length. The forward pass (preorder) passes the complementary sums down, so that every edge holds the sums of the tips on both of its sides, from which the root search evaluates any root position without another traversal.

v1: `fn clock_regression()` [packages/treetime/src/clock/clock_regression.rs#L161-L170](../../packages/treetime/src/clock/clock_regression.rs#L161-L170).
v0: [`packages/legacy/treetime/treetime/treeregression.py#L192-L283`](../../packages/legacy/treetime/treetime/treeregression.py#L192-L283).

Re-estimation mode: when `prev_clock_rate` is provided (during iterative timetree refinement), the regression uses solver-updated time lengths converted to divergence (`time_length * rate * gamma`) instead of input branch lengths (`fn edge_length()`, [packages/treetime/src/clock/clock_regression.rs#L315-L326](../../packages/treetime/src/clock/clock_regression.rs#L315-L326)).

v0 can weight the regression by branch variances (`--covariation`, `--tip-slack`); v1 removed these options ([kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md)).

Key functions: `clock_regression()`, `fit_clock_model()`.

---

## Exact Root Position

With the root at fraction $x$ of an edge of length $\ell$, the tips below the edge are $(1 - x)\ell$ and all others $x\ell$ further from the root. The date sums do not change and the divergence sums are linear or quadratic in $x$, so the score is a quadratic in $x$ and the parabola through $x = 0, \tfrac12, 1$ is exact; its vertex clamped to $[0, 1]$ is the optimum on that edge. On a terminal edge, a vertex past the tip gives $x = 1$, and the reroot then places a new root at the tip end instead of making the tip the root ([reroot.md](reroot.md), "Root Placement"). v0 brackets the split on a grid and refines it with bounded Brent minimization (`scipy.optimize.minimize_scalar`); v1 needs no 1D optimizer.

v1: `fn optimize_exact()` [packages/treetime/src/reroot/edge_search.rs#L93](../../packages/treetime/src/reroot/edge_search.rs#L93), cost function `struct BranchPointCostFunction` in the same file.

---

## IQD Outlier Detection

Identifies clock outliers using Tukey's fences (<a id="cite-5"></a>[Tukey 1977](https://doi.org/10.1007/978-1-4419-7976-6) [[5](#ref-5)]) applied to clock deviation residuals. The interquartile distance IQD = Q3 - Q1 measures the spread of the middle 50% of clock deviations. A leaf with `|clock_deviation| > threshold * IQD` is flagged as an outlier, where the threshold is controlled by `--clock-filter` (default 4.0, as in v0).

IQD-based detection is nonparametric and does not assume normally distributed residuals, making it appropriate for phylogenetic data where residual distributions can be skewed. For comparison, under a normal distribution IQD = 1.35 _ sigma, so `threshold _ IQD = 3 _ 1.35 _ sigma = 4.05 \* sigma` - roughly the 99.995% interval.

v1: [`packages/treetime/src/clock/clock_filter.rs#L22-L89`](../../packages/treetime/src/clock/clock_filter.rs#L22-L89).
v0: [`packages/legacy/treetime/treetime/clock_filter_methods.py#L5-L40`](../../packages/legacy/treetime/treetime/clock_filter_methods.py#L5-L40).

### References

- <a id="ref-5"></a>Tukey, John W. 1977. _Exploratory Data Analysis._ Addison-Wesley. ISBN 978-0-201-07616-5. [↩](#cite-5)

---

## Tree Rerooting

Edge inversion along the path from the old root to the new root position, with cleanup of trivial (single-child) nodes. Rerooting changes which node is the root of the tree without altering the unrooted topology. All edges along the old-root-to-new-root path are inverted (parent becomes child and vice versa), and the old root is removed if it becomes a trivial single-child node. An old root with one child would become a leaf after the inversion. When it has no date, the reroot removes it first and its child becomes the root, which is exact because a branch to a node without data contributes a factor of 1 to the likelihood and root-to-tip distances add along a path. A dated one becomes a dated leaf. `--keep-root` leaves the tree unchanged. A reroot that turns an internal node without a date into a leaf fails with an error ([packages/treetime/src/reroot/placement.rs](../../packages/treetime/src/reroot/placement.rs)). The clock fit (`fn fit_clock_model()`) only fits: every caller reroots first, then runs the backward regression pass on the new tree. The date-free methods (`min-dev`, `--outgroup`, `midpoint`) reroot once and keep the root for every later fit ([reroot.md](reroot.md), "Date-Based and Date-Free Methods").

v1: `fn reroot_at_best_root()` [packages/treetime/src/reroot/driver.rs#L16](../../packages/treetime/src/reroot/driver.rs#L16), topology move `fn apply_root_at_edge()` [packages/treetime/src/reroot/orchestrate.rs#L20](../../packages/treetime/src/reroot/orchestrate.rs#L20). See [reroot.md](reroot.md).

---

## Best Root Search

Every edge is searched for its exact optimum (above), and the minimum over the current root and all edges wins. With a required positive rate, an edge whose optimum has a non-positive rate is rejected; the rate is checked on the optimum rather than imposed as a constraint, so data without temporal signal stop with an error instead of fitting a rate near zero.

The positive-rate restriction is controlled by `force_positive_rate`. `--allow-negative-rate` relaxes it so negative-slope roots are accepted. When a non-positive rate survives to model construction (via `--keep-root`, `--allow-negative-rate`, or a root method without a slope search), the boundary policy differs by command: timetree never runs time inference on it (time inference needs `time = div / rate`) and takes the substitution rate of the Gaussian timetree instead, except in the pre-loop refits on a `least-squares` root, where it stops; the clock command warns and continues (it only reports the regression). See [decisions/timetree-runs-no-time-inference-on-a-non-positive-clock-rate.md](../decisions/timetree-runs-no-time-inference-on-a-non-positive-clock-rate.md).

The root search of the outlier pre-filter of the clock command always requires a positive rate, also with `--allow-negative-rate`, as in v0: `fn TreeTime.clock_filter()` reroots with the default `force_positive=True` ([packages/legacy/treetime/treetime/treetime.py#L486-L489](../../packages/legacy/treetime/treetime/treetime.py#L486-L489)), and `fn TreeRegression.find_best_root()` then accepts only roots with a non-negative slope ([packages/legacy/treetime/treetime/treeregression.py#L347-L353](../../packages/legacy/treetime/treetime/treeregression.py#L347-L353)). The flag applies to the root search after the filter only. When no root of the pre-filter has a positive rate, the error names the outlier filter and `--clock-filter=0` as the way to skip it.

v1: `fn estimate_clock_model_with_prefilter()` [packages/treetime/src/clock/pipeline.rs#L93-L136](../../packages/treetime/src/clock/pipeline.rs#L93-L136). The outlier filter of timetree has no root search of its own: it uses the clock fit of the first reroot, which follows `--allow-negative-rate` ([packages/treetime/src/timetree/pre_loop.rs#L45-L53](../../packages/treetime/src/timetree/pre_loop.rs#L45-L53)); see [kb/issues/M-timetree-allow-negative-rate-has-no-v0-counterpart.md](../issues/M-timetree-allow-negative-rate-has-no-v0-counterpart.md).

v1: `fn search_edges()` [packages/treetime/src/reroot/edge_search.rs#L19](../../packages/treetime/src/reroot/edge_search.rs#L19), positive-rate check `RootObjective::has_positive_rate()` [packages/treetime/src/reroot/params.rs](../../packages/treetime/src/reroot/params.rs).

---

## Additional Algorithms

- Root-to-Tip Distance: [`packages/treetime/src/clock/rtt.rs`](../../packages/treetime/src/clock/rtt.rs) - Tree traversal accumulating branch lengths from root to each leaf, producing the divergence values used in clock regression.
- Clock Model: [`packages/treetime/src/clock/clock_model.rs`](../../packages/treetime/src/clock/clock_model.rs) - Rate and intercept with `ClockModelStats::Estimated` (chi-squared, R-value and, when known, `rate_std`) or `ClockModelStats::Fixed` (from `--clock-rate`). A timetree rate from the marginal-likelihood step counts as estimated.
- Gaussian timetree: [`packages/treetime/src/clock/gaussian.rs`](../../packages/treetime/src/clock/gaussian.rs) - Gaussian message passing at a given rate: node-time marginals and the rate and root-time standard deviations.
- Date Constraints: [`packages/treetime/src/clock/date_constraints.rs`](../../packages/treetime/src/clock/date_constraints.rs) - `load_date_constraints()` (`#load_date_constraints`) assigns date distributions to tree nodes with validation.

---

## Unimplemented

See [unimplemented](unimplemented.md) for full details:

- Local outlier filter (z-score based)
- Full covariance matrix computation (TreeRegression), intentionally removed with `--covariation`
- Numerical Hessian for root position uncertainty

---

## File Index

| File                                                                                                       | Algorithms                                               |
| ---------------------------------------------------------------------------------------------------------- | -------------------------------------------------------- |
| [`packages/treetime/src/clock/clock_set.rs`](../../packages/treetime/src/clock/clock_set.rs)               | Least-squares sums                                       |
| [`packages/treetime/src/clock/clock_regression.rs`](../../packages/treetime/src/clock/clock_regression.rs) | Tree message passing, edge divergence, clock fit         |
| [`packages/treetime/src/clock/clock_model.rs`](../../packages/treetime/src/clock/clock_model.rs)           | Clock model and its statistics                           |
| [`packages/treetime/src/clock/gaussian.rs`](../../packages/treetime/src/clock/gaussian.rs)                 | Gaussian timetree, rate and root uncertainty             |
| [`packages/treetime/src/clock/clock_filter.rs`](../../packages/treetime/src/clock/clock_filter.rs)         | IQD outlier detection                                    |
| [`packages/treetime/src/reroot/`](../../packages/treetime/src/reroot/)                                     | Root search, exact optimum on an edge, tree rerooting    |
| [`packages/treetime/src/clock/date_constraints.rs`](../../packages/treetime/src/clock/date_constraints.rs) | Date constraint loading and validation                   |
| [`packages/treetime/src/clock/rtt.rs`](../../packages/treetime/src/clock/rtt.rs)                           | Root-to-tip distance computation                         |
| [`packages/treetime/src/clock/pipeline.rs`](../../packages/treetime/src/clock/pipeline.rs)                 | Clock pipeline                                           |
| [`packages/app-commands/src/commands/clock/run.rs`](../../packages/app-commands/src/commands/clock/run.rs) | Clock command orchestration                              |
