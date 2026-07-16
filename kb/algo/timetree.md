# Timetree Inference Algorithms

[Back to index](README.md)

## Belief Propagation

Two-pass message passing for time inference on the phylogenetic tree (<a id="cite-1"></a>[Pearl 1988](https://doi.org/10.1016/B978-0-08-051489-5.50001-5) [[1](#ref-1)]). The backward pass convolves and multiplies child distributions from leaves toward the root. The forward pass divides and convolves parent distributions from root toward leaves, refining each node's time estimate with information from the rest of the tree.

v1: [`packages/treetime/src/timetree/inference/backward_pass.rs`](../../packages/treetime/src/timetree/inference/backward_pass.rs), [`packages/treetime/src/timetree/inference/forward_pass.rs`](../../packages/treetime/src/timetree/inference/forward_pass.rs).
v0: [`packages/legacy/treetime/treetime/node_interpolator.py`](../../packages/legacy/treetime/treetime/node_interpolator.py).

- `propagate_distributions_backward()` (`#propagate_distributions_backward`) [packages/treetime/src/timetree/inference/backward_pass.rs#L17-L31](../../packages/treetime/src/timetree/inference/backward_pass.rs#L17-L31): skips bad branches (outlier and dateless leaves) so they do not constrain parent time
- `propagate_distributions_forward()` (`#propagate_distributions_forward`) [packages/treetime/src/timetree/inference/forward_pass.rs#L11-L22](../../packages/treetime/src/timetree/inference/forward_pass.rs#L11-L22): preserves internal node times when forward pass yields `None`

---

## ML Branch-Length Pre-Step

Before time inference, the timetree pipeline runs one pass of per-edge ML branch-length optimization to seed the belief propagation with better branch lengths than raw input values. This matches v0's `optimize_tree(max_iter=1)` calls (<a id="cite-3"></a>[Sagulenko, Puller, and Neher 2018](https://doi.org/10.1093/ve/vex042) [[3](#ref-3)], pipeline description).

v1: `optimize_branch_lengths_pre_step()` in [`packages/treetime/src/commands/timetree/run.rs`](../../packages/treetime/src/commands/timetree/run.rs), called twice: before and after rerooting.
v0: `optimize_tree(max_iter=1)` in [`packages/legacy/treetime/treetime/treetime.py`](../../packages/legacy/treetime/treetime/treetime.py) at lines 243 and 266.

The pre-step reuses `run_optimize_mixed()` from the optimize command infrastructure via trait object coercion (`PartitionTimetreeAll` extends `PartitionOptimizeOps`). Uses `BrentSqrt` method (matching v0's Brent in sqrt(t) space).

---

## Branch Length Distributions

Computes per-edge time distributions from partition contributions, clock rate, and gamma rate multiplier. This is a TreeTime-specific transformation that converts the branch-length likelihood (from sequence data) into a time-domain distribution suitable for belief propagation.

v1: [`packages/treetime/src/timetree/inference/branch_length_likelihood.rs`](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs).

- `compute_branch_length_distribution()` (`#compute_branch_length_distribution`) [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L44-L82](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L44-L82): evaluates the Poisson indel log-likelihood alongside the substitution log-likelihood on the grid (matching `run_optimize_mixed()`), then converts the branch-length grid to a time grid using `effective_clock_rate = clock_rate * gamma`, where `gamma > 1` means faster evolution (shorter time for same substitutions)

---

## Kingman Coalescent

The Kingman coalescent (<a id="cite-2"></a>[Kingman 1982](<https://doi.org/10.1016/0304-4149(82)90011-4>) [[2](#ref-2)]) provides a prior on node times. With $k(t)$ lineages and coalescent time scale $T_c(t)$, the total merger rate is $\lambda(t)=k(t)(k(t)-1)/(2T_c(t))$ and the per-lineage rate is $\kappa(t)=(k(t)-1)/(2T_c(t))$.

V1 expresses every quantity in decimal calendar years and stores the cumulative hazard $H(t)=\int_t^P\kappa(s)\,ds$, where $P$ is the most recent event. For a branch from parent $p$ to child $c$, the negative-log survival term is $H(t_p)-H(t_c)$. Grouping all branch terms by node gives:

| Piece                           | Neg-log value           | Purpose                                                                        |
| ------------------------------- | ----------------------- | ------------------------------------------------------------------------------ |
| Internal node (all, incl. root) | $H(t)-\ln\lambda(t)$     | Merger density and parent-side branch survival                                |
| Leaf                            | $-H(t_\mathrm{leaf})$     | Child-side survival credit                                                     |
| Root correction                 | $+H(t_\mathrm{root})$     | Root has no parent to supply a child-side subtraction                          |

These pieces sum to the edge-derived Kingman objective. An internal node with $m$ children receives multiplicity $m-1$, matching v0. Edge scoring distributes that merger-density term over the parent's $m$ outgoing edges, preserving the actual parent multiplicity correction.

v1: [`packages/treetime/src/coalescent/`](../../packages/treetime/src/coalescent/).
v0: [`packages/legacy/treetime/treetime/merger_models.py`](../../packages/legacy/treetime/treetime/merger_models.py).

The shared pipeline is:

- `collect_tree_events()` [packages/treetime/src/coalescent/events.rs](../../packages/treetime/src/coalescent/events.rs) collects calendar-dated sample and merger events.
- `compute_lineage_count_distribution()` [packages/treetime/src/coalescent/lineage_dynamics.rs](../../packages/treetime/src/coalescent/lineage_dynamics.rs) builds calendar-coordinate $k(t)$ with tested breakpoint sidedness.
- `compute_integral_merger_rate()` [packages/treetime/src/coalescent/integration.rs](../../packages/treetime/src/coalescent/integration.rs) computes $H(t)$ using the existing midpoint quadrature.
- `compute_coalescent_model()` [packages/treetime/src/coalescent/coalescent.rs](../../packages/treetime/src/coalescent/coalescent.rs) constructs one immutable model for inference and objective evaluation.

The backward pass combines every child message before evaluating one leaf, internal, or root cost on the completed distribution's existing coordinates. It subtracts the minimum finite cost before exponentiation and normalizes the result. Leaf factors affect only outgoing messages, so observed date distributions are not overwritten or compounded by repeated passes.

The first timetree pass runs without coalescent to establish node time distributions via backward+forward belief propagation. Coalescent contributions are computed from these established times on the second pass.
Detailed ownership and objective identities are documented in [kb/algo/coalescent-contribution-refactor.md](coalescent-contribution-refactor.md).

---

## Tc Optimization

Optimizes the coalescent time scale Tc in log space over the bracket [-20, 2] using Brent's method (<a id="cite-4"></a>[Brent 1973](https://doi.org/10.1007/978-3-0348-5952-3) [[4](#ref-4)]). Brent's method is a hybrid of parabolic interpolation and golden section search, achieving superlinear convergence without requiring derivatives.

`optimize_tc()` (`#optimize_tc`) [packages/treetime/src/coalescent/optimize_tc.rs](../../packages/treetime/src/coalescent/optimize_tc.rs) precomputes lineage state and inferred edge endpoints, then binds each constant candidate $T_c$ to `CoalescentModel` and minimizes the shared edge-derived objective. Nodes without inferred dates are skipped with a warning; reversed finite endpoints are errors.

Tc is re-optimized each iteration (from iteration 2 onward) using constant Tc. In skyline mode, constant Tc is used during loop iterations; full skyline fit is deferred to post-convergence.

---

## Skyline Coalescent

Piecewise-varying Tc(t) estimation with smoothness regularization, extending the constant-Tc coalescent to capture population size changes over time.

The classic skyline plot (<a id="cite-5"></a>[Pybus, Rambaut, and Harvey 2000](https://doi.org/10.1093/genetics/155.3.1429) [[5](#ref-5)]) estimates a step function for N_e(t) from inter-coalescent intervals. The MLE per interval is `N_hat_k = C(k) * w_k` where C(k) = k(k-1)/2 is the number of lineage pairs and w_k is the interval duration. This is noisy because each interval has at most one coalescent event.

The Bayesian skyline plot (<a id="cite-6"></a>[Drummond et al. 2005](https://doi.org/10.1093/molbev/msi103) [[6](#ref-6)]) reduces noise by grouping neighboring coalescent events into segments with shared N*e parameters, using an autocorrelated exponential prior linking adjacent population sizes. The skyride (<a id="cite-7"></a>[Minin, Bloomquist, and Suchard 2008](https://doi.org/10.1093/molbev/msn090) [[7](#ref-7)]) replaces grouping with a Gaussian Markov Random Field (GMRF) prior that penalizes changes per unit time: `pi(gamma|tau) ~ tau^((m-1)/2) * exp(-tau/2 _ sum((gamma_{k+1} - gamma_k)^2 / delta_k))` where tau is a global precision parameter and delta_k are time-aware weights.

TreeTime's skyline implementation uses a smoothness penalty similar to the GMRF but optimizes via Nelder-Mead rather than MCMC:

`optimize_skyline()` (`#optimize_skyline`) [packages/treetime/src/coalescent/skyline.rs#L68-L157](../../packages/treetime/src/coalescent/skyline.rs#L68-L157) minimizes:

```
cost = -log_likelihood + stiffness * sum(diff(log_Tc)^2) + regularization * boundary_penalty
```

The stiffness term penalizes rapid changes in log(Tc) between adjacent grid points (analogous to the GMRF precision). The boundary penalty discourages log(Tc) values outside [-100, 0].

`build_tc_distribution()` (`#build_tc_distribution`) [packages/treetime/src/coalescent/skyline.rs#L214-L256](../../packages/treetime/src/coalescent/skyline.rs#L214-L256) creates a `Distribution::Formula` with piecewise linear interpolation via binary search, producing a lazy Tc(t) function for the backward pass.

v1: [`packages/treetime/src/coalescent/skyline.rs`](../../packages/treetime/src/coalescent/skyline.rs).
v0: [`packages/legacy/treetime/treetime/merger_models.py#L281`](../../packages/legacy/treetime/treetime/merger_models.py#L281).
v1 uses Nelder-Mead (via `argmin` crate); v0 uses SLSQP (via `scipy.optimize.minimize`). See [known issue](../issues/M-timetree-skyline-nelder-mead-optimizer.md).

Skyline is re-optimized after the main iteration loop converges with stabilized node times, then a final timetree pass runs with the optimized Tc(t).

---

## Relaxed Clock

Relaxed molecular clocks allow the substitution rate to vary across branches, relaxing the strict-clock assumption (<a id="cite-9"></a>[Zuckerkandl and Pauling 1965](https://doi.org/10.1016/B978-1-4832-2734-4.50017-6) [[9](#ref-9)]) of a single uniform rate. Rate variation arises from differing generation times, population sizes, metabolic rates, and selective pressures across lineages.

TreeTime implements an autocorrelated model where descendant branch rates correlate with parent rates, following <a id="cite-10"></a>[Thorne, Kishino, and Painter 1998](https://doi.org/10.1093/oxfordjournals.molbev.a025892) [[10](#ref-10)]. Each branch carries a rate multiplier gamma: `effective_rate = clock_rate * gamma`. The model penalizes deviation from gamma=1 (slack parameter) and rate differences between parent and child (coupling parameter).

When coupling > 0, closely related branches have similar rates (autocorrelated clock, matching the biological expectation that rate-influencing traits evolve gradually). When coupling = 0, the model degenerates to an uncorrelated clock where each branch rate is independent (<a id="cite-11"></a>[Drummond et al. 2006](https://doi.org/10.1371/journal.pbio.0040088) [[11](#ref-11)]). <a id="cite-12"></a>[Lepage et al. 2007](https://doi.org/10.1093/molbev/msm193) [[12](#ref-12)] compared relaxed clock models and found that autocorrelated models perform better when rate variation is driven by lineage-specific traits, while uncorrelated models handle episodic rate changes better.

v1: [`packages/treetime/src/timetree/optimization/relaxed_clock.rs`](../../packages/treetime/src/timetree/optimization/relaxed_clock.rs).
v0: `TreeTime.relaxed_clock()` in [`packages/legacy/treetime/treetime/treetime.py`](../../packages/legacy/treetime/treetime/treetime.py).
CLI: `--relax <SLACK> <COUPLING>` (defaults 1.0, 1.0).

### Algorithm

`apply_relaxed_clock()` (`#apply_relaxed_clock`) [packages/treetime/src/timetree/optimization/relaxed_clock.rs#L25-L125](../../packages/treetime/src/timetree/optimization/relaxed_clock.rs#L25-L125) runs two passes:

- Postorder pass (lines 36-81): computes quadratic penalty coefficients k1, k2 per node. The penalty function is `stiffness * (gamma * actual_len - optimal_len)^2 + slack * (gamma - 1)^2`, with a coupling term `coupling * (gamma - gamma_child)^2`.
- Preorder pass (lines 86-114): computes optimal gamma per branch. Root: `gamma = max(0.1, -0.5 * k1 / k2)`. Non-root: `gamma = max(0.1, (coupling * parent_gamma - 0.5 * k1) / (coupling + k2))`.

`actual_len` and `optimal_len` are both substitutions per site. v1 converts its year-valued edge `time_length` with `clock_rate` before evaluating the penalty, matching v0's branch-length-valued `clock_length`.

The `one_mutation` parameter (sum of sequence lengths across all partitions) sets the scale for branch length penalties. Gamma values are clamped to a minimum of 0.1 to prevent degenerate solutions.

### Gamma consumption

`compute_branch_length_distribution()` (`#compute_branch_length_distribution`) uses `effective_clock_rate = clock_rate * gamma` to convert the branch-length grid to a time grid. `edge_divergence()` (`#edge_divergence`) uses `time_length * rate * gamma` when re-estimating divergence from solver-updated time lengths.

### References

- Zuckerkandl & Pauling (1962). "Molecular Disease, Evolution, and Genic Heterogeneity." In Kasha & Pullman (eds.), Horizons in Biochemistry, pp. 189-225. Academic Press

---

## Polytomy Resolution

A polytomy (multifurcation) is a node with more than two children, arising from insufficient phylogenetic signal to resolve the true bifurcating topology (soft polytomy) or from genuine simultaneous divergence (hard polytomy, rare). Tree builders (IQ-TREE, FastTree, RAxML) resolve zero-length branches into arbitrary bifurcations. TreeTime collapses these back into polytomies and re-resolves them in a way consistent with the temporal ordering of nodes.

v1 implements greedy deterministic resolution. v0 also supports stochastic coalescent-based resolution (not yet ported, see [unimplemented](unimplemented.md#stochastic-polytomy-resolution)).

v1: [`packages/treetime/src/timetree/optimization/polytomy.rs`](../../packages/treetime/src/timetree/optimization/polytomy.rs).

### Greedy algorithm

The algorithm iterates over all nodes with >2 children. For each polytomy, it computes pairwise likelihood gains for merging each pair of children under a new internal node, selects the pair with the highest gain, and merges them. This repeats until no pair exceeds the resolution threshold (default 0.05 log-likelihood units) or the node becomes binary.

This approach is deterministic and reproducible but biases toward caterpillar-like topologies: after the first merge creates a new internal node, subsequent merges preferentially attach to it (because it has the most informative branch distribution), creating an imbalanced subtree ([[3](#ref-3)], Section 2.6).

- `resolve_polytomies()` (`#resolve_polytomies`) [packages/treetime/src/timetree/optimization/polytomy.rs#L27-L32](../../packages/treetime/src/timetree/optimization/polytomy.rs#L27-L32): entry point with default threshold (0.05).
- `compute_merge_gain()` (`#compute_merge_gain`) [packages/treetime/src/timetree/optimization/polytomy.rs#L225](../../packages/treetime/src/timetree/optimization/polytomy.rs#L225): uses Brent optimization (via `argmin` crate) to find the optimal merge time and cost gain for a child pair.
- `merge_children()` (`#merge_children`) [packages/treetime/src/timetree/optimization/polytomy.rs#L345](../../packages/treetime/src/timetree/optimization/polytomy.rs#L345): creates a new internal node, adds parent-to-new-node edge, reparents the two children.
- `prepare_tree_after_topology_change()` (`#prepare_tree_after_topology_change`) [packages/treetime/src/timetree/optimization/polytomy.rs](../../packages/treetime/src/timetree/optimization/polytomy.rs): clears cached internal-node distributions and resets topology-dependent edge distributions, clock messages, and relaxed-clock rates; leaf date constraints, leaf `bad_branch` flags, branch lengths, and time lengths are preserved.

After resolution, partition data is reconciled via `reconcile_topology()` to add entries for new nodes/edges.

### Known issues

The zero-branch penalty for newly created internal branches differs from v0: v1 uses bare time difference, v0 scales by `gtr.mu * data.full_length`.

---

## Clock Filter / Outlier Detection

IQD-based outlier detection that marks leaves with anomalous root-to-tip divergence as bad branches.

v1: [`packages/treetime/src/timetree/optimization/clock_filter.rs`](../../packages/treetime/src/timetree/optimization/clock_filter.rs).

- `apply_outlier_bad_branches()` (`#apply_outlier_bad_branches`) [packages/treetime/src/timetree/optimization/clock_filter.rs#L75-L95](../../packages/treetime/src/timetree/optimization/clock_filter.rs#L75-L95): sets `bad_branch=true` on outlier leaves, then postorder propagation marks internal nodes bad only when all children are bad

---

## Convergence Monitoring

Tracks optimization loop convergence via sequence change counts and likelihood components. TreeTime-specific metrics.

v1: [`packages/treetime/src/timetree/convergence/`](../../packages/treetime/src/timetree/convergence/) (3 files).

- `TimetreeOptimizer` (`#TimetreeOptimizer`) [packages/treetime/src/timetree/convergence/metrics.rs#L16-L108](../../packages/treetime/src/timetree/convergence/metrics.rs#L16-L108): iteration controller that records `ConvergenceMetrics` per round and stops when converged or max iterations reached. Convergence criterion: `n_diff == 0 && n_resolved == 0`. Supports CSV tracelog output. Convergence can be suppressed (for skyline mode where constant Tc is used during iterations).
- `count_sequence_changes()` (`#count_sequence_changes`) [packages/treetime/src/timetree/convergence/sequence_changes.rs#L19-L47](../../packages/treetime/src/timetree/convergence/sequence_changes.rs#L19-L47): compares ancestral state snapshots position-by-position across all internal nodes and partitions
- `capture_ancestral_states()` (`#capture_ancestral_states`) [packages/treetime/src/timetree/convergence/sequence_changes.rs#L50-L70](../../packages/treetime/src/timetree/convergence/sequence_changes.rs#L50-L70): snapshots reconstructed sequences for all internal nodes. Captured before polytomy resolution to avoid inflating n_diff with newly created nodes.

### Likelihood components

- `compute_sequence_likelihood()` (`#compute_sequence_likelihood`) [packages/treetime/src/timetree/convergence/likelihood.rs#L11-L25](../../packages/treetime/src/timetree/convergence/likelihood.rs#L11-L25): sum of per-partition root log-likelihoods from marginal reconstruction
- `compute_positional_likelihood()` (`#compute_positional_likelihood`) [packages/treetime/src/timetree/convergence/likelihood.rs#L35-L76](../../packages/treetime/src/timetree/convergence/likelihood.rs#L35-L76): sum of log-probabilities of branch length distributions evaluated at inferred time durations. **v1-specific metric** - v0's `positional_LH` sums node-level marginal log-likelihoods from the forward pass. Both trend in the same direction during convergence but produce different numerical values.
- `compute_coalescent_likelihood()` (`#compute_coalescent_likelihood`) [packages/treetime/src/timetree/convergence/likelihood.rs#L80-L91](../../packages/treetime/src/timetree/convergence/likelihood.rs#L80-L91): total coalescent log-likelihood via `compute_coalescent_total_lh()`. Sums per-edge Kingman coalescent costs under the active Tc distribution.
- `compute_coalescent_total_lh()` (`#compute_coalescent_total_lh`) [packages/treetime/src/coalescent/total_lh.rs](../../packages/treetime/src/coalescent/total_lh.rs): constructs `CoalescentModel`, collects inferred endpoint data, and calls `sum_coalescent_cost()`.
- `collect_coalescent_edges()` (`#collect_coalescent_edges`) [packages/treetime/src/coalescent/edge_data.rs](../../packages/treetime/src/coalescent/edge_data.rs): collects inferred child and parent calendar dates plus actual parent multiplicity. Missing dates are skipped with a warning; reversed endpoints are rejected.
- `sum_coalescent_cost()` (`#sum_coalescent_cost`) [packages/treetime/src/coalescent/edge_data.rs](../../packages/treetime/src/coalescent/edge_data.rs): sums the model's endpoint-derived edge costs. The grouped merger-density share across all child edges equals the corresponding node contribution.

---

## Iterative EM-like Refinement

Alternates sequence reconstruction (E-step) and time inference (M-step), iterating until convergence ([[3](#ref-3)], Section 2.4). Each iteration optionally applies relaxed clock rate estimation, resolves polytomies, and re-estimates the clock model.

v1: [`packages/treetime/src/commands/timetree/refinement.rs`](../../packages/treetime/src/commands/timetree/refinement.rs).

- `run_refinement_iteration()` (`#run_refinement_iteration`) [packages/treetime/src/commands/timetree/refinement.rs#L21-L106](../../packages/treetime/src/commands/timetree/refinement.rs#L21-L106): per-iteration logic: relaxed clock, polytomy resolution, ancestral reconstruction, timetree inference, clock re-estimation. Captures ancestral state snapshots before polytomy resolution.

---

## Confidence Intervals

Node date uncertainty has two independent sources ([[3](#ref-3)], Section 2.2):

1. Mutation stochasticity. The Poisson process of substitution accumulation creates branch length uncertainty, which propagates through the backward/forward belief propagation passes. The marginal posterior distribution at each node captures this. Nodes constrained by many descendant dates have narrow posteriors; weakly constrained nodes have wide posteriors.
2. Clock rate uncertainty. The regression slope has a standard error from the 2x2 Hessian inverse (`ClockModel::cov()`). All node times scale inversely with the rate, so rate uncertainty propagates to all dates. Nodes near the root have the highest sensitivity: a 10% rate error shifts the root date by 10% of the tree depth.

v1: [`packages/treetime/src/timetree/confidence.rs`](../../packages/treetime/src/timetree/confidence.rs).
v0: `get_max_posterior_region()` (`#get_max_posterior_region`), `calc_rate_susceptibility()` (`#calc_rate_susceptibility`), `date_uncertainty_due_to_rate()` (`#date_uncertainty_due_to_rate`), `combine_confidence()` (`#combine_confidence`) in [`packages/legacy/treetime/treetime/clock_tree.py#L1010-L1230`](../../packages/legacy/treetime/treetime/clock_tree.py#L1010-L1230).

Both v0 and v1 use 90% confidence regions throughout: `CI_FRACTION=0.9` in v1, `fraction=0.9` in v0's `get_max_posterior_region`. The derived quantile bounds are `(0.05, 0.95)` (5% in each tail). The two sources use different statistical operations to produce "90% CI":

- Mutation: 90% HPD (highest posterior density) - the shortest interval containing 90% of the probability mass. Tighter than equal-tailed CI for skewed distributions.
- Rate: 90% equal-tailed interval via probit scaling - maps ±1 sigma rate sensitivity dates to the 5th/95th percentile using standard normal z-scores (z=±1.645).
- Combined: quadrature sum of deviations - a Gaussian-approximation heuristic, not a formal CI of the joint distribution.

### Orchestrator

`extract_confidence_intervals()` (`#extract_confidence_intervals`) computes per-node CI:

1. Mutation contribution from `Distribution::hpd_region(0.9)` on the marginal posterior
2. Rate contribution from `date_uncertainty_due_to_rate(dates, (0.05, 0.95))` on rate susceptibility triples
3. Combination via `combine_confidence()` quadrature sum, clipped to distribution domain
4. Delta distributions and nodes without data yield identity interval `[date, date]`

### Mutation contribution: HPD region

`Distribution::hpd_region(fraction)` in [packages/treetime-distribution/src/distribution_core/distribution.rs](../../packages/treetime-distribution/src/distribution_core/distribution.rs) finds the narrowest interval containing `fraction` of the probability mass:

- Interior peak: bisection on probability threshold. For each candidate threshold, the interval between left and right crossings of the PDF captures a certain probability mass. Bisection finds the threshold where the mass equals `fraction`. Assumes unimodality (returns a single contiguous interval).
- Boundary peak: one-sided quantile fallback. Left boundary: `[t_min, quantile(fraction)]`. Right boundary: `[quantile(1-fraction), t_max]`.
- Uniform/Range: centered interval of width `fraction * range`.
- Formula: equal-tailed fallback `confidence_interval(p_lo, 1-p_lo)`.

v0: `get_max_posterior_region(node, fraction=0.9)` (`clock_tree.py:1146-1230`) uses `scipy.optimize.minimize_scalar` (Brent) for the threshold search. v1 uses bisection, which is simpler and sufficient since the target function (CDF mass vs threshold) is monotonic.

### Rate contribution: susceptibility analysis

`compute_rate_susceptibility()` (`#compute_rate_susceptibility`) scales per-edge gamma values at `rate ± rate_std`, runs timetree inference three times, stores per-node date triples `[lower_date, central_date, upper_date]` sorted by date value.

v0: `calc_rate_susceptibility()` (`#calc_rate_susceptibility`) in [`clock_tree.py#L1010-L1066`](../../packages/legacy/treetime/treetime/clock_tree.py#L1010-L1066). Same gamma scaling approach: `gamma_new = gamma_orig * new_rate / current_rate`.

`date_uncertainty_due_to_rate()` (`#date_uncertainty_due_to_rate`) converts the date triple to CI bounds via the probit function (inverse normal CDF):

```
z_lower = probit(0.05) = -1.645
z_upper = probit(0.95) = +1.645
ci_lower = central + z_lower * |lower - central|
ci_upper = central + z_upper * |upper - central|
```

This is a split-normal (two-piece normal) approximation: lower and upper deviations are scaled independently as one-sided standard deviations. The standard normal probit is used for both sides. This matches v0 (`clock_tree.py:1083-1085`) and is valid when `rate_std << rate` (the rate-to-date mapping is approximately linear and the rate estimate is approximately Gaussian from regression CLT). For extreme asymmetry (rate_std approaching rate, triggering the `0.1 * rate` floor on the lower rate), the approximation is biased because the proper split-normal quantile function includes a normalization factor `A = sqrt(2/pi) / (sigma_l + sigma_u)` that the probit mapping omits.

v0 note: `date_uncertainty_due_to_rate` has default parameter `interval=(0.05, 0.095)` (`clock_tree.py:1068`), where `0.095` is a typo for `0.95`. The default is never reached because all callers pass explicit intervals.

### Combination

`combine_confidence()` (`#combine_confidence`) combines independent CI contributions via quadrature sum, treating the two deviations as independent Gaussian-like errors:

```
lower = center - sqrt((rate_lo - center)^2 + (mutation_lo - center)^2)
upper = center + sqrt((rate_hi - center)^2 + (mutation_hi - center)^2)
```

Clipped to distribution domain (physical limits). This is a first-order Gaussian approximation ([[3](#ref-3)], Section 2.2). The combined interval is wider than either source alone but narrower than the sum of deviations.

v0: `combine_confidence()` (`clock_tree.py:1090-1101`). Same formula.

### References

---

## Timetree Runner

Single timetree inference pass: branch distributions, coalescent model construction, backward pass, forward pass. This is the inner loop called by the estimation pipeline and by rate susceptibility analysis.

v1: [`packages/treetime/src/timetree/inference/runner.rs`](../../packages/treetime/src/timetree/inference/runner.rs).

- `run_timetree()` (`#run_timetree`) [packages/treetime/src/timetree/inference/runner.rs](../../packages/treetime/src/timetree/inference/runner.rs): computes branch distributions, constructs one optional calendar-coordinate coalescent model, then performs backward and forward propagation. Validates that the clock rate is positive.

---

## Estimation Pipeline

End-to-end timetree estimation from input to output.

v1: [`packages/treetime/src/commands/timetree/run.rs`](../../packages/treetime/src/commands/timetree/run.rs).

- `run_timetree_estimation()` (`#run_timetree_estimation`) [packages/treetime/src/commands/timetree/run.rs#L32-L253](../../packages/treetime/src/commands/timetree/run.rs#L32-L253): loads data, estimates clock model, reroots, runs clock filter, initializes partitions, runs initial timetree pass (without coalescent), optimizes Tc, runs second pass (with coalescent), iterates refinement with convergence monitoring, post-convergence skyline optimization, final marginal reconstruction for CI extraction, writes outputs.

---

## Reroot

Clock-based rerooting with partition state update. Searches for the root position that maximizes temporal signal (positive clock rate, minimal regression residuals).

v1: [`packages/treetime/src/timetree/optimization/reroot.rs`](../../packages/treetime/src/timetree/optimization/reroot.rs).

- `reroot_tree()` (`#reroot_tree`) [packages/treetime/src/timetree/optimization/reroot.rs#L20-L93](../../packages/treetime/src/timetree/optimization/reroot.rs#L20-L93): performs clock regression with rerooting, applies topology changes to partitions (edge split, edge merge, inverted edges), recomputes marginal messages

---

## References

- <a id="ref-1"></a>Pearl, Judea. 1988. _Probabilistic Reasoning in Intelligent Systems: Networks of Plausible Inference._ Morgan Kaufmann. ISBN 978-0-934613-73-2. [↩](#cite-1)
- <a id="ref-2"></a>Kingman, J. F. C. 1982. "The Coalescent." _Stochastic Processes and their Applications_ 13(3):235-248. https://doi.org/10.1016/0304-4149(82)90011-4 [↩](#cite-2)
- <a id="ref-3"></a>Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4(1):vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-3)
- <a id="ref-4"></a>Brent, Richard P. 1973. _Algorithms for Minimization Without Derivatives._ Prentice-Hall. ISBN 978-0-13-022335-7. [↩](#cite-4)
- <a id="ref-5"></a>Pybus, Oliver G., Andrew Rambaut, and Paul H. Harvey. 2000. "An Integrated Framework for the Inference of Viral Population History from Reconstructed Genealogies." _Genetics_ 155(3):1429-1437. https://doi.org/10.1093/genetics/155.3.1429 [↩](#cite-5)
- <a id="ref-6"></a>Drummond, Alexei J., Andrew Rambaut, Beth Shapiro, and Oliver G. Pybus. 2005. "Bayesian Coalescent Inference of Past Population Dynamics from Molecular Sequences." _Molecular Biology and Evolution_ 22(5):1185-1192. https://doi.org/10.1093/molbev/msi103 [↩](#cite-6)
- <a id="ref-7"></a>Minin, Vladimir N., Erik W. Bloomquist, and Marc A. Suchard. 2008. "Smooth Skyride Through a Rough Skyline: Bayesian Coalescent-Based Inference of Population Dynamics." _Molecular Biology and Evolution_ 25(7):1459-1471. https://doi.org/10.1093/molbev/msn090 [↩](#cite-7)
- <a id="ref-8"></a>Strimmer, Korbinian, and Oliver G. Pybus. 2001. "Exploring the Demographic History of DNA Sequences Using the Generalized Skyline Plot." _Molecular Biology and Evolution_ 18(12):2298-2305. https://doi.org/10.1093/oxfordjournals.molbev.a003776
- <a id="ref-9"></a>Zuckerkandl, Emile, and Linus Pauling. 1965. "Evolutionary Divergence and Convergence in Proteins." In _Evolving Genes and Proteins,_ edited by Vernon Bryson and Henry J. Vogel, 97-166. Academic Press. https://doi.org/10.1016/B978-1-4832-2734-4.50017-6 [↩](#cite-9)
- <a id="ref-10"></a>Thorne, Jeffrey L., Hirohisa Kishino, and Ian S. Painter. 1998. "Estimating the Rate of Evolution of the Rate of Molecular Evolution." _Molecular Biology and Evolution_ 15(12):1647-1657. https://doi.org/10.1093/oxfordjournals.molbev.a025892 [↩](#cite-10)
- <a id="ref-11"></a>Drummond, Alexei J., Simon Y. W. Ho, Matthew J. Phillips, and Andrew Rambaut. 2006. "Relaxed Phylogenetics and Dating with Confidence." _PLoS Biology_ 4(5):e88. https://doi.org/10.1371/journal.pbio.0040088 [↩](#cite-11)
- <a id="ref-12"></a>Lepage, Thomas, David Bryant, Herve Philippe, and Nicolas Lartillot. 2007. "A General Comparison of Relaxed Molecular Clock Models." _Molecular Biology and Evolution_ 24(12):2669-2680. https://doi.org/10.1093/molbev/msm193 [↩](#cite-12)
- <a id="ref-13"></a>Felsenstein, Joseph. 1985. "Confidence Limits on Phylogenies: An Approach Using the Bootstrap." _Evolution_ 39(4):783-791. https://doi.org/10.2307/2408678

---

## File Index

| File                                                                                                                   | Algorithms                                                            |
| ---------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| [`packages/treetime/src/timetree/inference/`](../../packages/treetime/src/timetree/inference/)       | Belief propagation, branch distributions, timetree runner             |
| [`packages/treetime/src/coalescent/`](../../packages/treetime/src/coalescent/)     | Kingman coalescent, skyline, Tc optimization                          |
| [`packages/treetime/src/timetree/optimization/`](../../packages/treetime/src/timetree/optimization/) | Polytomy, relaxed clock, reroot, clock filter                         |
| [`packages/treetime/src/timetree/convergence/`](../../packages/treetime/src/timetree/convergence/)   | Convergence monitoring, likelihood tracking, sequence change counting |
| [`packages/treetime/src/commands/timetree/output/`](../../packages/treetime/src/commands/timetree/output/)             | Confidence intervals, date output, plots                              |
| [`packages/treetime/src/commands/timetree/refinement.rs`](../../packages/treetime/src/commands/timetree/refinement.rs) | EM-like iterative refinement                                          |
| [`packages/treetime/src/commands/timetree/run.rs`](../../packages/treetime/src/commands/timetree/run.rs)               | End-to-end estimation pipeline                                        |
