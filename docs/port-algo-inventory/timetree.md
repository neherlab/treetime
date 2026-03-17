# Timetree Inference Algorithms

[Back to index](_index.md)

## Belief Propagation

Two-pass message passing for time inference on the phylogenetic tree (Pearl 1988). The backward pass convolves and multiplies child distributions from leaves toward the root. The forward pass divides and convolves parent distributions from root toward leaves, refining each node's time estimate with information from the rest of the tree.

v1: [`packages/treetime/src/commands/timetree/inference/backward_pass.rs`](../../packages/treetime/src/commands/timetree/inference/backward_pass.rs), [`packages/treetime/src/commands/timetree/inference/forward_pass.rs`](../../packages/treetime/src/commands/timetree/inference/forward_pass.rs).
v0: [`packages/legacy/treetime/treetime/node_interpolator.py`](../../packages/legacy/treetime/treetime/node_interpolator.py).

- `propagate_distributions_backward()` (`#propagate_distributions_backward`) [packages/treetime/src/commands/timetree/inference/backward_pass.rs#L17-L31](../../packages/treetime/src/commands/timetree/inference/backward_pass.rs#L17-L31): skips bad branches (outlier and dateless leaves) so they do not constrain parent time
- `propagate_distributions_forward()` (`#propagate_distributions_forward`) [packages/treetime/src/commands/timetree/inference/forward_pass.rs#L11-L22](../../packages/treetime/src/commands/timetree/inference/forward_pass.rs#L11-L22): preserves internal node times when forward pass yields `None`

Reference: Pearl (1988). "Probabilistic Reasoning in Intelligent Systems." Morgan Kaufmann.

---

## Branch Length Distributions

Computes per-edge time distributions from partition contributions, clock rate, and gamma rate multiplier. This is a TreeTime-specific transformation that converts the branch-length likelihood (from sequence data) into a time-domain distribution suitable for belief propagation.

v1: [`packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs`](../../packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs).

- `compute_branch_length_distribution()` (`#compute_branch_length_distribution`) [packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L31-L63](../../packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L31-L63): converts branch length grid to time grid using `effective_clock_rate = clock_rate * gamma`, where `gamma > 1` means faster evolution (shorter time for same substitutions)

---

## Kingman Coalescent

The Kingman coalescent (Kingman 1982) provides a prior on internal node times based on population genetics. Under neutral evolution, k lineages in a population of effective size N_e merge backward in time at a pairwise rate, producing a total coalescence rate `lambda(t) = k(k-1) / (2*Tc(t))` where Tc is the coalescence time scale (proportional to N_e). The per-lineage merger rate is `kappa(t) = (k(t)-1) / (2*Tc(t))`, and the integral merger rate `I(t) = integral kappa(t') dt'` accumulates the probability of no coalescence up to time t.

At an internal node with m children, the coalescent contribution to the node's time distribution has two parts: a survival factor `exp(-I(t))` for leaves (probability of no coalescence), and a merger density `lambda(t)^(m-1) * exp(-I(t)*(m-1))` for internal nodes. In neg-log space: leaf contribution = `-I(t)`, internal = `multiplicity * (I(t) - log(lambda(t)))`.

v1: [`packages/treetime/src/commands/timetree/coalescent/`](../../packages/treetime/src/commands/timetree/coalescent/) (9 files).
v0: [`packages/legacy/treetime/treetime/merger_models.py`](../../packages/legacy/treetime/treetime/merger_models.py).

The coalescent contribution pipeline chains four steps:

- `collect_tree_events()` (`#collect_tree_events`) [packages/treetime/src/commands/timetree/coalescent/events.rs#L13-L53](../../packages/treetime/src/commands/timetree/coalescent/events.rs#L13-L53): traverses the graph breadth-first, collecting `(time, delta_branches)` tuples. Leaves get +1, internal nodes with k children get -(k-1).
- `compute_lineage_count_distribution()` (`#compute_lineage_count_distribution`) [packages/treetime/src/commands/timetree/coalescent/lineage_dynamics.rs#L14-L41](../../packages/treetime/src/commands/timetree/coalescent/lineage_dynamics.rs#L14-L41): aggregates events into a `PiecewiseConstantFn` representing k(t).
- `compute_integral_merger_rate()` (`#compute_integral_merger_rate`) [packages/treetime/src/commands/timetree/coalescent/integration.rs#L44-L94](../../packages/treetime/src/commands/timetree/coalescent/integration.rs#L44-L94): computes I(t) = integral of kappa(t). Clamping `k_clamped = max(0.5, k-1)` matches v0.
- `compute_node_contributions()` (`#compute_node_contributions`) [packages/treetime/src/commands/timetree/coalescent/contributions.rs#L26-L71](../../packages/treetime/src/commands/timetree/coalescent/contributions.rs#L26-L71): leaf contribution = -I(t), internal = multiplicity \* (I(t) - log(lambda(t))).

`compute_coalescent_contributions()` (`#compute_coalescent_contributions`) [packages/treetime/src/commands/timetree/coalescent/coalescent.rs#L61-L79](../../packages/treetime/src/commands/timetree/coalescent/coalescent.rs#L61-L79) orchestrates the full pipeline.

The first timetree pass runs without coalescent to establish node time distributions via backward+forward belief propagation. Coalescent contributions are computed from these established times on the second pass.

References:

- Kingman (1982). "The coalescent." Stochastic Processes and Applications, 13(3):235-248. doi:10.1016/0304-4149(82)90011-4
- Sagulenko, Puller & Neher (2018). "TreeTime." Virus Evolution, 4(1):vex042. doi:10.1093/ve/vex042

---

## Tc Optimization

Optimizes the coalescent time scale Tc in log space over the bracket [-20, 2] using Brent's method (Brent 1973). Brent's method is a hybrid of parabolic interpolation and golden section search, achieving superlinear convergence without requiring derivatives.

`optimize_tc()` (`#optimize_tc`) [packages/treetime/src/commands/timetree/coalescent/optimize_tc.rs#L44-L84](../../packages/treetime/src/commands/timetree/coalescent/optimize_tc.rs#L44-L84) precomputes lineage counts and per-node branch data, then minimizes negative total coalescent likelihood. The cost function (`TcCostFunction`) builds a constant `Distribution` at each evaluation, computes `compute_integral_merger_rate()`, and sums per-branch costs. Nodes with undetermined branch length are skipped with a warning.

Tc is re-optimized each iteration (from iteration 2 onward) using constant Tc. In skyline mode, constant Tc is used during loop iterations; full skyline fit is deferred to post-convergence.

Reference: Brent (1973). "Algorithms for Minimization Without Derivatives." Prentice-Hall, Chapter 5.

---

## Skyline Coalescent

Piecewise-varying Tc(t) estimation with smoothness regularization, extending the constant-Tc coalescent to capture population size changes over time.

The classic skyline plot (Pybus, Rambaut & Harvey 2000) estimates a step function for N_e(t) from inter-coalescent intervals. The MLE per interval is `N_hat_k = C(k) * w_k` where C(k) = k(k-1)/2 is the number of lineage pairs and w_k is the interval duration. This is noisy because each interval has at most one coalescent event.

The Bayesian skyline plot (Drummond, Rambaut, Shapiro & Pybus 2005) reduces noise by grouping neighboring coalescent events into segments with shared N*e parameters, using an autocorrelated exponential prior linking adjacent population sizes. The skyride (Minin, Bloomquist & Suchard 2008) replaces grouping with a Gaussian Markov Random Field (GMRF) prior that penalizes changes per unit time: `pi(gamma|tau) ~ tau^((m-1)/2) * exp(-tau/2 _ sum((gamma_{k+1} - gamma_k)^2 / delta_k))` where tau is a global precision parameter and delta_k are time-aware weights.

TreeTime's skyline implementation uses a smoothness penalty similar to the GMRF but optimizes via Nelder-Mead rather than MCMC:

`optimize_skyline()` (`#optimize_skyline`) [packages/treetime/src/commands/timetree/coalescent/skyline.rs#L68-L157](../../packages/treetime/src/commands/timetree/coalescent/skyline.rs#L68-L157) minimizes:

```
cost = -log_likelihood + stiffness * sum(diff(log_Tc)^2) + regularization * boundary_penalty
```

The stiffness term penalizes rapid changes in log(Tc) between adjacent grid points (analogous to the GMRF precision). The boundary penalty discourages log(Tc) values outside [-100, 0].

`build_tc_distribution()` (`#build_tc_distribution`) [packages/treetime/src/commands/timetree/coalescent/skyline.rs#L214-L256](../../packages/treetime/src/commands/timetree/coalescent/skyline.rs#L214-L256) creates a `Distribution::Formula` with piecewise linear interpolation via binary search, producing a lazy Tc(t) function for the backward pass.

v1: [`packages/treetime/src/commands/timetree/coalescent/skyline.rs`](../../packages/treetime/src/commands/timetree/coalescent/skyline.rs).
v0: [`packages/legacy/treetime/treetime/merger_models.py#L281`](../../packages/legacy/treetime/treetime/merger_models.py#L281).
v1 uses Nelder-Mead (via `argmin` crate); v0 uses SLSQP (via `scipy.optimize.minimize`). See [known issue](../port-known-issues/M-timetree-skyline-nelder-mead-optimizer.md).

Skyline is re-optimized after the main iteration loop converges with stabilized node times, then a final timetree pass runs with the optimized Tc(t).

References:

- Pybus, Rambaut & Harvey (2000). "An integrated framework for the inference of viral population history from reconstructed genealogies." Genetics, 155(3):1429-1437.
- Drummond, Rambaut, Shapiro & Pybus (2005). "Bayesian coalescent inference of past population dynamics from molecular sequences." Mol Biol Evol, 22(5):1185-1192. doi:10.1093/molbev/msi103
- Minin, Bloomquist & Suchard (2008). "Smooth skyride through a rough skyline: Bayesian coalescent-based inference of population dynamics." Mol Biol Evol, 25(7):1459-1471. doi:10.1093/molbev/msn090
- Strimmer & Pybus (2001). "Exploring the demographic history of DNA sequences using the generalized skyline plot." Mol Biol Evol, 18(12):2298-2305.

---

## Relaxed Clock

Relaxed molecular clocks allow the substitution rate to vary across branches, relaxing the strict-clock assumption (Zuckerkandl & Pauling 1962) of a single uniform rate. Rate variation arises from differing generation times, population sizes, metabolic rates, and selective pressures across lineages.

TreeTime implements an autocorrelated model where descendant branch rates correlate with parent rates, following Thorne, Kishino & Painter (1998). Each branch carries a rate multiplier gamma: `effective_rate = clock_rate * gamma`. The model penalizes deviation from gamma=1 (slack parameter) and rate differences between parent and child (coupling parameter).

When coupling > 0, closely related branches have similar rates (autocorrelated clock, matching the biological expectation that rate-influencing traits evolve gradually). When coupling = 0, the model degenerates to an uncorrelated clock where each branch rate is independent (Drummond et al. 2006). Lepage et al. (2007) compared relaxed clock models and found that autocorrelated models perform better when rate variation is driven by lineage-specific traits, while uncorrelated models handle episodic rate changes better.

v1: [`packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs`](../../packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs).
v0: `TreeTime.relaxed_clock()` in [`packages/legacy/treetime/treetime/treetime.py`](../../packages/legacy/treetime/treetime/treetime.py).
CLI: `--relax <SLACK> <COUPLING>` (defaults 1.0, 1.0).

### Algorithm

`apply_relaxed_clock()` (`#apply_relaxed_clock`) [packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs#L25-L125](../../packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs#L25-L125) runs two passes:

- **Postorder pass** (lines 36-81): computes quadratic penalty coefficients k1, k2 per node. The penalty function is `stiffness * (gamma * actual_len - optimal_len)^2 + slack * (gamma - 1)^2`, with a coupling term `coupling * (gamma - gamma_child)^2`.
- **Preorder pass** (lines 86-114): computes optimal gamma per branch. Root: `gamma = max(0.1, -0.5 * k1 / k2)`. Non-root: `gamma = max(0.1, (coupling * parent_gamma - 0.5 * k1) / (coupling + k2))`.

The `one_mutation` parameter (sum of sequence lengths across all partitions) sets the scale for branch length penalties. Gamma values are clamped to a minimum of 0.1 to prevent degenerate solutions.

### Gamma consumption

`compute_branch_length_distribution()` (`#compute_branch_length_distribution`) uses `effective_clock_rate = clock_rate * gamma` to convert the branch-length grid to a time grid. `edge_divergence()` (`#edge_divergence`) uses `time_length * rate * gamma` when re-estimating divergence from solver-updated time lengths.

### References

- Thorne, Kishino & Painter (1998). "Estimating the rate of evolution of the rate of molecular evolution." Mol Biol Evol, 15(12):1647-1657. doi:10.1093/oxfordjournals.molbev.a025892
- Drummond, Ho, Phillips & Rambaut (2006). "Relaxed phylogenetics and dating with confidence." PLOS Biology, 4(5):e88. doi:10.1371/journal.pbio.0040088
- Lepage, Bryant, Philippe & Lartillot (2007). "A general comparison of relaxed molecular clock models." Mol Biol Evol, 24(12):2669-2680. doi:10.1093/molbev/msm193

---

## Polytomy Resolution

A polytomy (multifurcation) is a node with more than two children, arising from insufficient phylogenetic signal to resolve the true bifurcating topology (soft polytomy) or from genuine simultaneous divergence (hard polytomy, rare). Tree builders (IQ-TREE, FastTree, RAxML) resolve zero-length branches into arbitrary bifurcations. TreeTime collapses these back into polytomies and re-resolves them in a way consistent with the temporal ordering of nodes.

v1 implements greedy deterministic resolution. v0 also supports stochastic coalescent-based resolution (not yet ported, see [unimplemented](unimplemented.md#stochastic-polytomy-resolution)).

v1: [`packages/treetime/src/commands/timetree/optimization/polytomy.rs`](../../packages/treetime/src/commands/timetree/optimization/polytomy.rs).

### Greedy algorithm

The algorithm iterates over all nodes with >2 children. For each polytomy, it computes pairwise likelihood gains for merging each pair of children under a new internal node, selects the pair with the highest gain, and merges them. This repeats until no pair exceeds the resolution threshold (default 0.05 log-likelihood units) or the node becomes binary.

This approach is deterministic and reproducible but biases toward caterpillar-like topologies: after the first merge creates a new internal node, subsequent merges preferentially attach to it (because it has the most informative branch distribution), creating an imbalanced subtree (Sagulenko et al. 2018, Section 2.6).

- `resolve_polytomies()` (`#resolve_polytomies`) [packages/treetime/src/commands/timetree/optimization/polytomy.rs#L27-L32](../../packages/treetime/src/commands/timetree/optimization/polytomy.rs#L27-L32): entry point with default threshold (0.05).
- `compute_merge_gain()` (`#compute_merge_gain`) [packages/treetime/src/commands/timetree/optimization/polytomy.rs#L225](../../packages/treetime/src/commands/timetree/optimization/polytomy.rs#L225): uses Brent optimization (via `argmin` crate) to find the optimal merge time and cost gain for a child pair.
- `merge_children()` (`#merge_children`) [packages/treetime/src/commands/timetree/optimization/polytomy.rs#L345](../../packages/treetime/src/commands/timetree/optimization/polytomy.rs#L345): creates a new internal node, adds parent-to-new-node edge, reparents the two children.
- `prepare_tree_after_topology_change()` (`#prepare_tree_after_topology_change`) [packages/treetime/src/commands/timetree/optimization/polytomy.rs#L432-L454](../../packages/treetime/src/commands/timetree/optimization/polytomy.rs#L432-L454): clears cached distributions on internal nodes; leaf date constraints and `bad_branch` flags are preserved.

After resolution, partition data is reconciled via `reconcile_topology()` to add entries for new nodes/edges.

### Known issues

The zero-branch penalty for newly created internal branches differs from v0: v1 uses bare time difference, v0 scales by `gtr.mu * data.full_length`. See [known issue](../port-known-issues/M-timetree-polytomy-zero-branch-penalty.md).

---

## Clock Filter / Outlier Detection

IQD-based outlier detection that marks leaves with anomalous root-to-tip divergence as bad branches.

v1: [`packages/treetime/src/commands/timetree/optimization/clock_filter.rs`](../../packages/treetime/src/commands/timetree/optimization/clock_filter.rs).

- `apply_outlier_bad_branches()` (`#apply_outlier_bad_branches`) [packages/treetime/src/commands/timetree/optimization/clock_filter.rs#L75-L95](../../packages/treetime/src/commands/timetree/optimization/clock_filter.rs#L75-L95): sets `bad_branch=true` on outlier leaves, then postorder propagation marks internal nodes bad only when all children are bad

---

## Convergence Monitoring

Tracks optimization loop convergence via sequence change counts and likelihood components. TreeTime-specific metrics.

v1: [`packages/treetime/src/commands/timetree/convergence/`](../../packages/treetime/src/commands/timetree/convergence/) (3 files).

- `TimetreeOptimizer` (`#TimetreeOptimizer`) [packages/treetime/src/commands/timetree/convergence/metrics.rs#L16-L108](../../packages/treetime/src/commands/timetree/convergence/metrics.rs#L16-L108): iteration controller that records `ConvergenceMetrics` per round and stops when converged or max iterations reached. Convergence criterion: `n_diff == 0 && n_resolved == 0`. Supports CSV tracelog output. Convergence can be suppressed (for skyline mode where constant Tc is used during iterations).
- `count_sequence_changes()` (`#count_sequence_changes`) [packages/treetime/src/commands/timetree/convergence/sequence_changes.rs#L19-L47](../../packages/treetime/src/commands/timetree/convergence/sequence_changes.rs#L19-L47): compares ancestral state snapshots position-by-position across all internal nodes and partitions
- `capture_ancestral_states()` (`#capture_ancestral_states`) [packages/treetime/src/commands/timetree/convergence/sequence_changes.rs#L50-L70](../../packages/treetime/src/commands/timetree/convergence/sequence_changes.rs#L50-L70): snapshots reconstructed sequences for all internal nodes. Captured before polytomy resolution to avoid inflating n_diff with newly created nodes.

### Likelihood components

- `compute_sequence_likelihood()` (`#compute_sequence_likelihood`) [packages/treetime/src/commands/timetree/convergence/likelihood.rs#L11-L25](../../packages/treetime/src/commands/timetree/convergence/likelihood.rs#L11-L25): sum of per-partition root log-likelihoods from marginal reconstruction
- `compute_positional_likelihood()` (`#compute_positional_likelihood`) [packages/treetime/src/commands/timetree/convergence/likelihood.rs#L35-L76](../../packages/treetime/src/commands/timetree/convergence/likelihood.rs#L35-L76): sum of log-probabilities of branch length distributions evaluated at inferred time durations. **v1-specific metric** - v0's `positional_LH` sums node-level marginal log-likelihoods from the forward pass. Both trend in the same direction during convergence but produce different numerical values.
- `compute_coalescent_likelihood()` (`#compute_coalescent_likelihood`) [packages/treetime/src/commands/timetree/convergence/likelihood.rs#L81-L83](../../packages/treetime/src/commands/timetree/convergence/likelihood.rs#L81-L83): stub, returns `None`

---

## Iterative EM-like Refinement

Alternates sequence reconstruction (E-step) and time inference (M-step), iterating until convergence (Sagulenko et al. 2018, Section 2.4). Each iteration optionally applies relaxed clock rate estimation, resolves polytomies, and re-estimates the clock model.

v1: [`packages/treetime/src/commands/timetree/refinement.rs`](../../packages/treetime/src/commands/timetree/refinement.rs).

- `run_refinement_iteration()` (`#run_refinement_iteration`) [packages/treetime/src/commands/timetree/refinement.rs#L21-L106](../../packages/treetime/src/commands/timetree/refinement.rs#L21-L106): per-iteration logic: relaxed clock, polytomy resolution, ancestral reconstruction, timetree inference, clock re-estimation. Captures ancestral state snapshots before polytomy resolution.

Reference: Sagulenko, Puller & Neher (2018). "TreeTime." Virus Evolution, 4(1):vex042, Section 2.4.

---

## Confidence Intervals

Node date uncertainty has two independent sources (Sagulenko, Puller & Neher 2018, Section 2.5):

1. **Mutation stochasticity.** The Poisson process of substitution accumulation creates branch length uncertainty, which propagates through the backward/forward belief propagation passes. The marginal posterior distribution at each node captures this. Nodes constrained by many descendant dates have narrow posteriors; weakly constrained nodes have wide posteriors.
2. **Clock rate uncertainty.** The regression slope has a standard error from the 2x2 Hessian inverse (`ClockModel::cov()`). All node times scale inversely with the rate, so rate uncertainty propagates to all dates. Nodes near the root have the highest sensitivity: a 10% rate error shifts the root date by 10% of the tree depth.

v1 implements both sources. Source (2) is activated via `--vary-rate` or `--confidence` with `--covariation`/`--clock-std-dev`.

v1: [`packages/treetime/src/commands/timetree/output/confidence.rs`](../../packages/treetime/src/commands/timetree/output/confidence.rs).
v0: `get_confidence_interval()` (`#get_confidence_interval`), `get_max_posterior_region()` (`#get_max_posterior_region`), `calc_rate_susceptibility()` (`#calc_rate_susceptibility`) in [`packages/legacy/treetime/treetime/clock_tree.py#L1010-L1230`](../../packages/legacy/treetime/treetime/clock_tree.py#L1010-L1230).

### Mutation-derived CI (implemented)

`extract_confidence_intervals()` (`#extract_confidence_intervals`) [packages/treetime/src/commands/timetree/output/confidence.rs#L46-L71](../../packages/treetime/src/commands/timetree/output/confidence.rs#L46-L71) computes 95% CI for each node using `Distribution::quantile()` at 0.025 and 0.975. Delta distributions yield identity interval `[date, date]`.

v0 uses a different method: `get_max_posterior_region(n, fraction=0.9)` finds the narrowest 90% highest posterior density (HPD) region around the peak via `scipy.optimize.minimize_scalar`. HPD regions are narrower than equal-tailed quantile intervals for skewed distributions. The porting decision should clarify which approach to match.

### CI combination

`combine_confidence()` (`#combine_confidence`) in [packages/treetime/src/commands/timetree/output/confidence.rs](../../packages/treetime/src/commands/timetree/output/confidence.rs) combines independent CI contributions via quadrature sum:

```
lower = center - hypot(c1_lower - center, c2_lower - center)
upper = center + hypot(c1_upper - center, c2_upper - center)
```

Clipped to physical limits. Called by `extract_confidence_intervals()` when rate variation data is present.

### Rate susceptibility (implemented)

`compute_rate_susceptibility()` (`#compute_rate_susceptibility`) in [packages/treetime/src/commands/timetree/output/confidence.rs](../../packages/treetime/src/commands/timetree/output/confidence.rs). Scales per-edge gamma values at rate +/- 1 sigma (matching v0 gamma scaling approach), runs timetree 3 times, stores per-node date triples sorted by date. `date_uncertainty_due_to_rate()` (`#date_uncertainty_due_to_rate`) converts triples to CI via probit function (erfinv).

v0 reference: `calc_rate_susceptibility()` (`#calc_rate_susceptibility`) in [packages/legacy/treetime/treetime/clock_tree.py#L1010-L1066](../../packages/legacy/treetime/treetime/clock_tree.py#L1010-L1066).

#### CI method

Both v1 and v0 use 90% HPD (highest posterior density) regions for the mutation contribution. HPD finds the narrowest interval containing 90% of the probability mass, which is tighter than equal-tailed quantile intervals for skewed distributions. v1 falls back to one-sided quantile CI when the peak is at a boundary, matching v0's `get_max_posterior_region` boundary handling (`clock_tree.py:1175-1178`).

### References

- Sagulenko, Puller & Neher (2018). "TreeTime." Virus Evolution, 4(1):vex042, Section 2.5. doi:10.1093/ve/vex042
- Felsenstein (1985). "Confidence limits on phylogenies: an approach using the bootstrap." Evolution, 39(4):783-791. doi:10.2307/2408678

---

## Timetree Runner

Single timetree inference pass: branch distributions, coalescent contributions, backward pass, forward pass. This is the inner loop called by the estimation pipeline and by rate susceptibility analysis.

v1: [`packages/treetime/src/commands/timetree/inference/runner.rs`](../../packages/treetime/src/commands/timetree/inference/runner.rs).

- `run_timetree()` (`#run_timetree`) [packages/treetime/src/commands/timetree/inference/runner.rs#L23-L74](../../packages/treetime/src/commands/timetree/inference/runner.rs#L23-L74): computes branch distributions (from partitions or input lengths), optional coalescent contributions, then backward+forward propagation. Validates clock rate is positive.

---

## Estimation Pipeline

End-to-end timetree estimation from input to output.

v1: [`packages/treetime/src/commands/timetree/run.rs`](../../packages/treetime/src/commands/timetree/run.rs).

- `run_timetree_estimation()` (`#run_timetree_estimation`) [packages/treetime/src/commands/timetree/run.rs#L32-L253](../../packages/treetime/src/commands/timetree/run.rs#L32-L253): loads data, estimates clock model, reroots, runs clock filter, initializes partitions, runs initial timetree pass (without coalescent), optimizes Tc, runs second pass (with coalescent), iterates refinement with convergence monitoring, post-convergence skyline optimization, final marginal reconstruction for CI extraction, writes outputs.

---

## Reroot

Clock-based rerooting with partition state update. Searches for the root position that maximizes temporal signal (positive clock rate, minimal regression residuals).

v1: [`packages/treetime/src/commands/timetree/optimization/reroot.rs`](../../packages/treetime/src/commands/timetree/optimization/reroot.rs).

- `reroot_tree()` (`#reroot_tree`) [packages/treetime/src/commands/timetree/optimization/reroot.rs#L20-L93](../../packages/treetime/src/commands/timetree/optimization/reroot.rs#L20-L93): performs clock regression with rerooting, applies topology changes to partitions (edge split, edge merge, inverted edges), recomputes marginal messages

---

## File Index

| File                                                                                                                   | Algorithms                                                            |
| ---------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| [`packages/treetime/src/commands/timetree/inference/`](../../packages/treetime/src/commands/timetree/inference/)       | Belief propagation, branch distributions, timetree runner             |
| [`packages/treetime/src/commands/timetree/coalescent/`](../../packages/treetime/src/commands/timetree/coalescent/)     | Kingman coalescent, skyline, Tc optimization                          |
| [`packages/treetime/src/commands/timetree/optimization/`](../../packages/treetime/src/commands/timetree/optimization/) | Polytomy, relaxed clock, reroot, clock filter                         |
| [`packages/treetime/src/commands/timetree/convergence/`](../../packages/treetime/src/commands/timetree/convergence/)   | Convergence monitoring, likelihood tracking, sequence change counting |
| [`packages/treetime/src/commands/timetree/output/`](../../packages/treetime/src/commands/timetree/output/)             | Confidence intervals, date output, plots                              |
| [`packages/treetime/src/commands/timetree/refinement.rs`](../../packages/treetime/src/commands/timetree/refinement.rs) | EM-like iterative refinement                                          |
| [`packages/treetime/src/commands/timetree/run.rs`](../../packages/treetime/src/commands/timetree/run.rs)               | End-to-end estimation pipeline                                        |
