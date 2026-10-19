# Timetree Inference Algorithms

## Belief Propagation

Two-pass message passing for time inference on the phylogenetic tree (<a id="cite-1"></a>[Pearl 1988](https://doi.org/10.1016/C2009-0-27609-4) [[1](#ref-1)]). The backward pass convolves and multiplies child distributions from leaves toward the root. The forward pass divides and convolves parent distributions from root toward leaves, refining each node's time estimate with information from the rest of the tree.

v1: [`packages/treetime/src/timetree/inference/backward_pass.rs`](../../packages/treetime/src/timetree/inference/backward_pass.rs), [`packages/treetime/src/timetree/inference/forward_pass.rs`](../../packages/treetime/src/timetree/inference/forward_pass.rs).
v0: [`packages/legacy/treetime/treetime/node_interpolator.py`](../../packages/legacy/treetime/treetime/node_interpolator.py).

- `fn run_timetree` [packages/treetime/src/timetree/inference/runner.rs#L30](../../packages/treetime/src/timetree/inference/runner.rs#L30): runs one time inference. Its inputs are `struct TimeInferenceInputs` [packages/treetime/src/timetree/inference/runner.rs#L81-L91](../../packages/treetime/src/timetree/inference/runner.rs#L81-L91) and an optional coalescent prior
- `fn derive_bad_branches` call [packages/treetime/src/timetree/inference/runner.rs#L54](../../packages/treetime/src/timetree/inference/runner.rs#L54): derives the bad-branch flags from the date constraints and the leaf flags
- Branch likelihoods [packages/treetime/src/timetree/inference/runner.rs#L56-L65](../../packages/treetime/src/timetree/inference/runner.rs#L56-L65): one per edge, from the input branch lengths or from the marginal reconstruction
- Passes [packages/treetime/src/timetree/inference/runner.rs#L67-L71](../../packages/treetime/src/timetree/inference/runner.rs#L67-L71): the backward pass runs first, then the forward pass. The backward output stays local to the call
- `struct TimeInference` [packages/treetime/src/timetree/inference/result.rs#L12-L16](../../packages/treetime/src/timetree/inference/result.rs#L12-L16): the result, with one map per role: bad-branch flags, branch likelihoods, and node posteriors
- `struct NodePosterior` [packages/treetime/src/timetree/inference/result.rs#L58-L61](../../packages/treetime/src/timetree/inference/result.rs#L58-L61): the distribution, the committed time and the contradicted flag of one node
- Every input of a call comes from its arguments. No value from a previous call enters the next call
- `propagate_distributions_backward()` (`#propagate_distributions_backward`) [packages/treetime/src/timetree/inference/backward_pass.rs](../../packages/treetime/src/timetree/inference/backward_pass.rs): a bad node (outlier or dateless leaf, or an undated node whose children are all bad) sends no message, so it does not constrain its parent's time; multiplies each node's fixed input date constraint into the combined child messages. A node without child messages and without a date has no subtree distribution
- `propagate_distributions_forward()` (`#propagate_distributions_forward`) [packages/treetime/src/timetree/inference/forward_pass.rs](../../packages/treetime/src/timetree/inference/forward_pass.rs): the posterior of a node without subtree distribution is the message from its parent alone, as in v0 (`clock_tree.py#L890-L891`); otherwise the parent posterior divided by this node's message, convolved across the edge and multiplied by the subtree distribution

What separates a leaf from an internal node here is not that it is a leaf but whether its date is
exact. A node given an exact date is pinned to it: nothing to refine, and no projection onto the
parent time. Every other node -- internal, or a leaf whose date is uncertain, ranged, or missing --
has its time inferred, and the forward pass refines it against the message coming down from the
parent. The date read from the input is an input of every backward pass, and the posterior is a
separate output, so a refined posterior never feeds back as evidence. See
[kb/decisions/timetree-uncertain-leaf-dates-are-inferred.md](../decisions/timetree-uncertain-leaf-dates-are-inferred.md).

---

## ML Branch-Length Pre-Step

Before time inference, the timetree pipeline runs one pass of per-edge ML branch-length optimization to seed the belief propagation with better branch lengths than raw input values. This matches v0's `optimize_tree(max_iter=1)` calls (<a id="cite-3"></a>[Sagulenko, Puller, and Neher 2018](https://doi.org/10.1093/ve/vex042) [[3](#ref-3)], pipeline description).

v1: `fn optimize_branch_lengths` [packages/treetime/src/timetree/pre_loop.rs#L209](../../packages/treetime/src/timetree/pre_loop.rs#L209), called through `fn ml_optimize` [packages/treetime/src/timetree/pre_loop.rs#L115](../../packages/treetime/src/timetree/pre_loop.rs#L115) twice: before and after rerooting ([packages/treetime/src/timetree/pre_loop.rs#L27-L71](../../packages/treetime/src/timetree/pre_loop.rs#L27-L71)).
v0: `optimize_tree(max_iter=1)` in [`packages/legacy/treetime/treetime/treetime.py`](../../packages/legacy/treetime/treetime/treetime.py) at lines 243 and 266.

The pre-step gathers the per-edge contributions of the marginal reconstruction and calls `fn optimize_edge_lengths`, with a zero indel rate when indels are off [packages/treetime/src/timetree/pre_loop.rs#L220-L240](../../packages/treetime/src/timetree/pre_loop.rs#L220-L240). It uses the `Brent` method (matching v0's Brent in sqrt(t) space).

---

## Branch Length Distributions

Computes per-edge time distributions from partition contributions, clock rate, and gamma rate multiplier. This is a TreeTime-specific transformation that converts the branch-length likelihood (from sequence data) into a time-domain distribution suitable for belief propagation.

v1: [`packages/treetime/src/timetree/inference/branch_length_likelihood.rs`](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs).

- `compute_branch_length_distribution()` (`#compute_branch_length_distribution`) [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L17-L65](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L17-L65): evaluates the Poisson indel log-likelihood alongside the substitution log-likelihood on the grid (matching `optimize_edge_lengths()`), then converts the branch-length grid to a time grid using the argument `branch_rate` (the clock rate times `gamma`), where `gamma > 1` means faster evolution (shorter time for same substitutions)

---

## Kingman Coalescent

The Kingman coalescent (<a id="cite-2"></a>[Kingman 1982](<https://doi.org/10.1016/0304-4149(82)90011-4>) [[2](#ref-2)]) provides a prior on node times. With $k(t)$ lineages and coalescent time scale $T_c(t)$, the total merger rate is $\lambda(t)=k(t)(k(t)-1)/(2T_c(t))$ and the per-lineage rate is $\kappa(t)=(k(t)-1)/(2T_c(t))$.

V1 expresses every quantity in decimal calendar years and works with the expected lineage count $\bar k(t)=\sum_u w_u S_u(t)$ instead of the lineage count at point times: $w_u = 1 - n_{\text{children}}(u)$ and $S_u(t)=P(t_u>t)$ from the node posteriors (`LineageExpectation`). The prior is the mean-field linearization of the Kingman density, local to each node: $-\log f_v(x)=\int^x [w_v(w_v-1)/2+w_v\bar k_{-v}(t)]/T_c(t)\,dt+m_v\ln T_c(x)$ with $\bar k_{-v}=\bar k-w_vS_v$ and $m_v$ the mergers at $v$ (`CoalescentPrior::node_contribution`). The same expectation gives the $T_c$ fit, the polytomy merger rates, the reported coalescent log-likelihood and the convergence trace.

v1: [`packages/treetime/src/coalescent/mean_field.rs`](../../packages/treetime/src/coalescent/mean_field.rs).
v0: [`packages/legacy/treetime/treetime/merger_models.py`](../../packages/legacy/treetime/treetime/merger_models.py).

The backward pass combines every child message before applying the node's prior term on the completed distribution's existing coordinates, in negative-log space. A leaf applies its factor $f_v$ to its own subtree distribution, which is both the evidence the forward pass combines into the tip posterior and the message the leaf sends to its parent, so the factor enters the tip posterior once and the parent once. Each pass starts again from the date constraints, so repeated passes do not compound it. A tip whose date is a range gets a grid over the range before the factor is applied, so the factor is evaluated across the range and not only at its ends.

The init phase seeds $\bar k(t)$ from the Gaussian node-time marginals at the regression rate, infers the node times once under that prior and rebuilds $\bar k(t)$ from their posteriors; it is then held fixed for the rest of the run, because deriving the prior from the times it is used to infer makes the loop chase a receding target. $T_c$ and the reported coalescent likelihood keep reading the current posteriors. See [kb/decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md).

Detailed ownership and objective identities are documented in [kb/algo/coalescent-contribution-refactor.md](coalescent-contribution-refactor.md).

---

## Tc Optimization

Optimizes the coalescent time scale Tc in log space over the bracket [-20, 2] using Brent's method (<a id="cite-4"></a>[Brent 1973](https://maths-people.anu.edu.au/brent/pub/pub011.html) [[4](#ref-4)]). Brent's method is a hybrid of parabolic interpolation and golden section search, achieving superlinear convergence without requiring derivatives.

`estimate_coalescent_tc()` (`#estimate_coalescent_tc`) [packages/treetime/src/timetree/coalescent_timescale.rs#L40](../../packages/treetime/src/timetree/coalescent_timescale.rs#L40) estimates a constant $T_c$ as the one-segment skyline: it calls `optimize_skyline()` with `n_points = 1`, whose optimum is the closed form $T_c^{*}=I/M$ from the expected pair-rate integral $I$ and expected merger count $M$ of the posteriors (`LineageExpectation::segment_statistics`).

Tc is re-optimized each iteration (from iteration 2 onward) using constant Tc. In skyline mode, constant Tc is used during loop iterations; full skyline fit is deferred to post-convergence.

---

## Skyline Coalescent

Piecewise-varying Tc(t) estimation with smoothness regularization, extending the constant-Tc coalescent to capture population size changes over time.

The classic skyline plot (<a id="cite-5"></a>[Pybus, Rambaut, and Harvey 2000](https://doi.org/10.1093/genetics/155.3.1429) [[5](#ref-5)]) estimates a step function for N_e(t) from inter-coalescent intervals. The MLE per interval is `N_hat_k = C(k) * w_k` where C(k) = k(k-1)/2 is the number of lineage pairs and w_k is the interval duration. This is noisy because each interval has at most one coalescent event.

The Bayesian skyline plot (<a id="cite-6"></a>[Drummond et al. 2005](https://doi.org/10.1093/molbev/msi103) [[6](#ref-6)]) reduces noise by grouping neighboring coalescent events into segments with shared N*e parameters, using an autocorrelated exponential prior linking adjacent population sizes. The skyride (<a id="cite-7"></a>[Minin, Bloomquist, and Suchard 2008](https://doi.org/10.1093/molbev/msn090) [[7](#ref-7)]) replaces grouping with a Gaussian Markov Random Field (GMRF) prior that penalizes changes per unit time: `pi(gamma|tau) ~ tau^((m-1)/2) * exp(-tau/2 _ sum((gamma_{k+1} - gamma_k)^2 / delta_k))` where tau is a global precision parameter and delta_k are time-aware weights.

TreeTime's skyline implementation uses a smoothness penalty similar to the GMRF but finds the optimum by Newton steps rather than MCMC:

`optimize_skyline()` (`#optimize_skyline`) [packages/treetime/src/coalescent/skyline.rs#L15](../../packages/treetime/src/coalescent/skyline.rs#L15) minimizes, over $z_i=\ln T_{c,i}$ on equal segments of the span of the expected lineage counts,

```
cost = sum_i (I_i exp(-z_i) + M_i z_i) + stiffness / 2 * sum_i (z_{i+1} - z_i)^2
```

with $I_i$, $M_i$ the expected pair-rate integral and merger count of segment $i$ (`LineageExpectation::segment_statistics`). The data term is the negative expected complete-data log-likelihood; the stiffness term penalizes rapid changes in log(Tc) between adjacent segments (analogous to the GMRF precision).

`build_tc_distribution()` (`#build_tc_distribution`) [packages/treetime/src/coalescent/skyline.rs](../../packages/treetime/src/coalescent/skyline.rs) creates a piecewise-constant `Distribution::Formula` for coalescent scoring. `SkylineResult::tc_schedule` preserves the same step boundaries for consumers that integrate hazards.

v1: [`packages/treetime/src/coalescent/skyline.rs`](../../packages/treetime/src/coalescent/skyline.rs).
v0: [`packages/legacy/treetime/treetime/merger_models.py#L281`](../../packages/legacy/treetime/treetime/merger_models.py#L281).
v1 uses a convex Newton solve in log $T_c$; v0 uses SLSQP via `scipy.optimize.minimize`. See [kb/decisions/coalescent-skyline-convex-log-tc.md](../decisions/coalescent-skyline-convex-log-tc.md).

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

`apply_relaxed_clock()` (`#apply_relaxed_clock`) [packages/treetime/src/timetree/optimization/relaxed_clock.rs#L10-L96](../../packages/treetime/src/timetree/optimization/relaxed_clock.rs#L10-L96) runs two passes:

- Postorder pass ([packages/treetime/src/timetree/optimization/relaxed_clock.rs#L24-L58](../../packages/treetime/src/timetree/optimization/relaxed_clock.rs#L24-L58)): computes quadratic penalty coefficients k1, k2 per node. The penalty function is `stiffness * (gamma * actual_len - optimal_len)^2 + slack * (gamma - 1)^2`, with a coupling term `coupling * (gamma - gamma_child)^2`.
- Preorder pass ([packages/treetime/src/timetree/optimization/relaxed_clock.rs#L62-L89](../../packages/treetime/src/timetree/optimization/relaxed_clock.rs#L62-L89)): computes optimal gamma per branch. Root: `gamma = max(0.1, -0.5 * k1 / k2)`. Non-root: `gamma = max(0.1, (coupling * parent_gamma - 0.5 * k1) / (coupling + k2))`.

`actual_len` and `optimal_len` are both substitutions per site. v1 converts its year-valued edge `time_length` with `clock_rate` before evaluating the penalty. That is not the contrast v0 fits: v0 compares `clock_length` against `mutation_length`, whereas both v1 inputs are ML-derived. See [kb/issues/M-timetree-consumers-read-unconstrained-branch-lengths.md](../issues/M-timetree-consumers-read-unconstrained-branch-lengths.md).

The `one_mutation` parameter (sum of sequence lengths across all partitions) sets the scale for branch length penalties. Gamma values are clamped to a minimum of 0.1 to prevent degenerate solutions.

### Gamma consumption

`compute_branch_length_distribution()` (`#compute_branch_length_distribution`) uses the argument `branch_rate` (the clock rate times `gamma`, computed by the caller) to convert the branch-length grid to a time grid. `edge_divergence()` (`#edge_divergence`) uses `time_length * rate * gamma` when re-estimating divergence from solver-updated time lengths.

### References

- Zuckerkandl & Pauling (1962). "Molecular Disease, Evolution, and Genic Heterogeneity." In Kasha & Pullman (eds.), Horizons in Biochemistry, pp. 189-225. Academic Press

---

## Polytomy Resolution

A polytomy (multifurcation) is a node with more than two children, arising from insufficient phylogenetic signal to resolve the true bifurcating topology (soft polytomy) or from genuine simultaneous divergence (hard polytomy, rare). Tree builders (IQ-TREE, FastTree, RAxML) resolve zero-length branches into arbitrary bifurcations. v0 collapses these back into polytomies before the timetree loop (`prune_short_branches()`) and re-resolves them in a way consistent with the temporal ordering of nodes. v1 does not collapse them, so it re-resolves only polytomies already present in the input tree: [kb/issues/M-timetree-pre-loop-short-branch-pruning-missing.md](../issues/M-timetree-pre-loop-short-branch-pruning-missing.md).

v1 uses mutation-conditioned stochastic coalescent sampling. v0 provides stochastic and greedy generators, but v1 corrects three documented v0 defects in the stochastic rate selection, parent bound, and arrival handling.

v1: [`packages/treetime/src/timetree/optimization/polytomy/`](../../packages/treetime/src/timetree/optimization/polytomy/).
v0: `generate_subtree()` in [`packages/legacy/treetime/treetime/treetime.py`](../../packages/legacy/treetime/treetime/treetime.py).

For live lineages, let $M$ be the remaining substitution count and let $R$ be the mutation-free lineages. The competing event rates are

$$R_{\mathrm{mut}}=\mu M, \qquad R_{\mathrm{coal}}=\max(0, |R|-1)\kappa(t),$$

where $\mu$ is the whole-alignment mutation rate and $\kappa(t)=\max(1/2,\bar k(t)-1)/(2T_c(t))$ is the per-branch merger-rate schedule of the frozen expected lineage counts (`LineageExpectation::merger_rate_schedule`). One unit-exponential hazard threshold is integrated across lineage arrivals, merger-rate breakpoints, and the parent bound. A mutation event selects a lineage in exact proportion to its integer remaining count. A coalescent event selects two mutation-free lineages uniformly.

- `simulate_subtree()` in [packages/treetime/src/timetree/optimization/polytomy/sweep.rs](../../packages/treetime/src/timetree/optimization/polytomy/sweep.rs) performs pure seeded simulation and returns a validated merger plan. It rejects non-finite times and rates, rejects rate overflow, and leaves residual polytomies when the parent bound stops the sweep.
- `apply_plan()` in [packages/treetime/src/timetree/optimization/polytomy/apply.rs](../../packages/treetime/src/timetree/optimization/polytomy/apply.rs) validates the complete plan before graph mutation, then creates internal nodes and reparents children.
- `resolve_polytomies()` in [packages/treetime/src/timetree/optimization/polytomy/resolve.rs](../../packages/treetime/src/timetree/optimization/polytomy/resolve.rs) collects exact reconstructed substitution counts when available, applies the plan, removes obsolete single-child nodes, and rebuilds the graph.
- `resolve_polytomies()` returns the times of the nodes it creates. `require_internal_node_times()` checks that every internal node has a finite time, before the first graph mutation and again for the resolved tree.

After resolution, partition data is reconciled with the new topology before inference resumes. The same seed and input state produce the same sampled plan. See [kb/decisions/timetree-stochastic-polytomy-resolution.md](../decisions/timetree-stochastic-polytomy-resolution.md) and the three `timetree-stochastic-resolve-*` entries in [kb/v0-errata/](../v0-errata/).

---

## Clock Filter / Outlier Detection

IQD-based outlier detection that marks leaves with anomalous root-to-tip divergence as bad branches.

v1: the IQD filter in [`packages/treetime/src/clock/clock_filter.rs`](../../packages/treetime/src/clock/clock_filter.rs), the timetree outlier report and leaf flags in [`packages/treetime/src/timetree/optimization/outliers.rs`](../../packages/treetime/src/timetree/optimization/outliers.rs).

- `mark_outlier_leaves()` (`#mark_outlier_leaves`) [packages/treetime/src/timetree/optimization/outliers.rs](../../packages/treetime/src/timetree/optimization/outliers.rs): adds the outlier leaves to the leaf bad-branch flags, which start as the leaves without a date
- `derive_bad_branches()` (`#derive_bad_branches`) [packages/treetime/src/timetree/inference/bad_branches.rs](../../packages/treetime/src/timetree/inference/bad_branches.rs): derives the flags of the current tree at the start of every time inference. A node with its own date is never bad; any other internal node is bad when all its children are bad. See [kb/decisions/timetree-dated-internal-node-never-bad-branch.md](../decisions/timetree-dated-internal-node-never-bad-branch.md)

---

## Convergence Monitoring

Tracks optimization loop convergence via sequence change counts and likelihood components. TreeTime-specific metrics.

v1: [`packages/treetime/src/timetree/convergence/`](../../packages/treetime/src/timetree/convergence/).

- `struct TimetreeOptimizer` [packages/treetime/src/timetree/convergence/optimizer.rs#L16-L116](../../packages/treetime/src/timetree/convergence/optimizer.rs#L16-L116): iteration controller that records `ConvergenceMetrics` per round and stops when converged or max iterations reached. Each record also goes to an optional `TraceSink` (the CSV tracelog). Convergence can be suppressed (for skyline mode where constant Tc is used during iterations)
- `fn ConvergenceMetrics::has_converged` [packages/treetime/src/timetree/convergence/metrics.rs#L20-L23](../../packages/treetime/src/timetree/convergence/metrics.rs#L20-L23): a round converges when no polytomy was resolved (`n_resolved == 0`) and the node times settled: the largest node-time change is below `NODE_TIME_TOLERANCE_YEARS`, or `n_diff == 0` when no node time is available
- `count_sequence_changes()` (`#count_sequence_changes`) [packages/treetime/src/timetree/convergence/sequence_changes.rs#L8-L23](../../packages/treetime/src/timetree/convergence/sequence_changes.rs#L8-L23): compares ancestral state snapshots position-by-position across all internal nodes and partitions
- `capture_ancestral_states()` (`#capture_ancestral_states`) [packages/treetime/src/timetree/convergence/sequence_changes.rs#L25-L34](../../packages/treetime/src/timetree/convergence/sequence_changes.rs#L25-L34): snapshots reconstructed sequences for all internal nodes. Captured before polytomy resolution to avoid inflating n_diff with newly created nodes.

### Log-likelihood components

Components use `LogLh`, while distribution grids and probability arrays remain native `f64`. The tracelog CSV columns are the fields of `ConvergenceMetrics`: `log_lh_seq`, `log_lh_time` and `log_lh_coal`.

- `compute_sequence_log_lh()` (`#compute_sequence_log_lh`) [packages/treetime/src/timetree/convergence/likelihood.rs#L13-L24](../../packages/treetime/src/timetree/convergence/likelihood.rs#L13-L24): sum of per-partition root log-likelihoods from marginal reconstruction
- `log_lh_time`: `TimeBackward::log_marginal_likelihood()` (`#log_marginal_likelihood`) [packages/treetime/src/timetree/inference/result.rs#L18](../../packages/treetime/src/timetree/inference/result.rs#L18), stored in `TimeInference` by `run_timetree` [packages/treetime/src/timetree/inference/runner.rs#L103](../../packages/treetime/src/timetree/inference/runner.rs#L103): log mass of the root's subtree distribution, whose amplitude carries the levels of all messages and factors (see [distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md#amplitude-re-peaking-and-normalization)). Branch-length likelihoods enter normalized to their peaks, so the value is the log marginal likelihood relative to each branch's sequence optimum, coalescent prior included.
- `log_lh_coal`: `CoalescentTimescale::expected_log_likelihood()` (`#expected_log_likelihood`) [packages/treetime/src/timetree/coalescent_timescale.rs#L142](../../packages/treetime/src/timetree/coalescent_timescale.rs#L142): expected complete-data coalescent log-likelihood $-\sum_i I_i/T_{c,i}-\sum_i M_i\ln T_{c,i}$ of the round's posteriors under the round's $T_c$, the objective the $T_c$ fit maximizes

---

## Iterative EM-like Refinement

Alternates sequence reconstruction (E-step) and time inference (M-step), iterating until convergence ([[3](#ref-3)], Section 2.4). Each iteration optionally applies relaxed clock rate estimation, resolves polytomies, and re-estimates the clock model.

A single backward and forward pass over per-branch length likelihoods already gives a timetree. Iteration is needed only for the optional parts. Some form a cycle: node times determine the coalescent likelihood, and the coalescent prior changes the node times. Others alternate in EM style: branch times are optimized for a fixed coalescent model or fixed rate multipliers, then those are optimized for fixed times. v1 runs the loop up to `--max-iter` times even when no optional part is enabled.

v1: [`packages/treetime/src/timetree/round.rs`](../../packages/treetime/src/timetree/round.rs), driven by [`packages/treetime/src/timetree/refinement_loop.rs`](../../packages/treetime/src/timetree/refinement_loop.rs).

- `refinement_round()` (`#refinement_round`) [packages/treetime/src/timetree/round.rs](../../packages/treetime/src/timetree/round.rs): takes one round state by value and names each state transition: relaxed-clock update (`relax_clock`), topology refinement (`refine_topology`), inference rebuild (`refresh_times`), ancestral-state comparison, and clock re-estimation (`update_clock_model`; with a coalescent prior the rate steps toward the maximum of `log_lh_time`, see [timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md)). `TopologyOutcome` distinguishes an unchanged tree from a changed tree with a resolved-node count, so callers handle topology status explicitly. A changed topology triggers partition reconciliation and a reset of the relaxed-clock rate multipliers to 1 in `fn refine_topology` [packages/treetime/src/timetree/round.rs#L316-L317](../../packages/treetime/src/timetree/round.rs#L316-L317). `fn refresh_times` [packages/treetime/src/timetree/round.rs#L340-L359](../../packages/treetime/src/timetree/round.rs#L340-L359) then runs marginal reconstruction and time inference with the active coalescent prior. The round returns the latest `TimeInference` and rate multipliers as values.

---

## Confidence Intervals

Node date uncertainty has two independent sources ([[3](#ref-3)], Section 2.2):

1. Mutation stochasticity. The Poisson process of substitution accumulation creates branch length uncertainty, which propagates through the backward/forward belief propagation passes. The marginal posterior distribution at each node captures this. Nodes constrained by many descendant dates have narrow posteriors; weakly constrained nodes have wide posteriors.
2. Clock rate uncertainty. The rate standard deviation is `--clock-std-dev` when given, otherwise, for an estimated rate with a known sequence length, the one of the Gaussian timetree at the final rate (`fn gaussian_clock()`, [`packages/treetime/src/clock/gaussian.rs`](../../packages/treetime/src/clock/gaussian.rs); [kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md)). v0 takes it from the covariance of the `--covariation` regression. All node times scale inversely with the rate, so rate uncertainty propagates to all dates. Nodes near the root have the highest sensitivity: a 10% rate error shifts the root date by 10% of the tree depth.

v1: [`packages/treetime/src/timetree/confidence.rs`](../../packages/treetime/src/timetree/confidence.rs), [`packages/treetime/src/timetree/max_posterior_region.rs`](../../packages/treetime/src/timetree/max_posterior_region.rs).
v0: `get_max_posterior_region()` (`#get_max_posterior_region`), `calc_rate_susceptibility()` (`#calc_rate_susceptibility`), `date_uncertainty_due_to_rate()` (`#date_uncertainty_due_to_rate`), `combine_confidence()` (`#combine_confidence`) in [`packages/legacy/treetime/treetime/clock_tree.py#L1010-L1230`](../../packages/legacy/treetime/treetime/clock_tree.py#L1010-L1230).

Both v0 and v1 use 90% confidence regions throughout: `CI_FRACTION=0.9` in v1, `fraction=0.9` in v0's `get_max_posterior_region`. The derived quantile bounds are `(0.05, 0.95)` (5% in each tail). The two sources use different statistical operations to produce "90% CI":

- Mutation: 90% HPD (highest posterior density) - the shortest interval containing 90% of the probability mass. Tighter than equal-tailed CI for skewed distributions.
- Rate: 90% equal-tailed interval via probit scaling - maps ±1 sigma rate sensitivity dates to the 5th/95th percentile using standard normal z-scores (z=±1.645).
- Combined: quadrature sum of deviations - a Gaussian-approximation heuristic, not a formal CI of the joint distribution.

### Orchestrator

`extract_confidence_intervals()` (`#extract_confidence_intervals`) computes per-node CI:

1. Mutation contribution from `max_posterior_region()` (`#max_posterior_region`) on the marginal time posterior
2. Rate contribution from `date_uncertainty_due_to_rate()` on rate susceptibility triples, with the quantile pair that matches the peak location: `(0.05, 0.95)` for an interior peak, `(0, 0.9)` for a peak at the latest date, `(0.1, 1)` for a peak at the earliest date, as v0 passes them through `get_confidence_interval`
3. Combination via `confidence_bounds()` and `combine_in_quadrature()`: quadrature sum, clipped to distribution domain
4. Delta distributions yield `[date, date]`; nodes without a distribution and without rate dates yield the identity interval `[date, date]`

### Mutation contribution: HPD region

v0: `get_max_posterior_region(node, fraction=0.9)` (`clock_tree.py:1146-1230`) finds the narrowest interval containing `fraction` of the probability mass, using `scipy.optimize.minimize_scalar` (Brent) for the threshold search.

v1: `max_posterior_region()` (`#max_posterior_region`) converts the negative-log ordinates to peak-relative densities $e^{-(y - y_{\min})}$ and integrates them on the posterior grid with the trapezoid rule, as v0 does on its grid. The CDF is linear between grid points. The peak is the grid minimum of $y$; a tie resolves to the latest date, which is v0's first minimum on its time-before-present axis.

- Peak inside the grid: the region is bounded by the crossings of a level $h$ above the peak on each side, and $h$ is chosen so that the region holds 90% of the mass. The mass is piecewise linear in $h$, so v1 solves for $h$ exactly where v0 runs a Brent minimizer on $(0.9 - \text{mass}(h))^2$. On a shallow posterior that minimizer stops on a plateau and v0 reports the whole support ([kb/v0-errata/timetree-max-posterior-region-brent-plateau.md](../v0-errata/timetree-max-posterior-region-brent-plateau.md)). On a side that is not monotone, v1 uses the first crossing outward from the peak ([kb/v0-errata/timetree-max-posterior-region-multimodal-side.md](../v0-errata/timetree-max-posterior-region-multimodal-side.md))
- Peak at the latest date: `[Q(0.1), t_max]`; peak at the earliest date: `[t_min, Q(0.9)]`, where `Q` is the inverse CDF. v0 reaches these through `get_confidence_interval(node, (0, 0.9))` and `(0.1, 1)` on its reversed axis
- Point distribution: `[t, t]`

Tail mass beyond the grid is not integrated, as in v0. On `ebola/362` and `rsv/a/500` the posterior tails beyond the grid hold at most `2.3e-5` of the mass, and almost all posteriors have hard boundaries on both sides.

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

`combine_in_quadrature()` (`#combine_in_quadrature`) combines independent CI contributions via quadrature sum, treating the two deviations as independent Gaussian-like errors:

```
lower = center - sqrt((rate_lo - center)^2 + (mutation_lo - center)^2)
upper = center + sqrt((rate_hi - center)^2 + (mutation_hi - center)^2)
```

Clipped to distribution domain (physical limits). This is a first-order Gaussian approximation ([[3](#ref-3)], Section 2.2). The combined interval is wider than either source alone but narrower than the sum of deviations.

v0: `combine_confidence()` (`clock_tree.py:1090-1101`). Same formula.

### References

---

## Timetree Runner

Single timetree inference pass: branch distributions, backward pass, forward pass. This is the inner loop called by the estimation pipeline and by rate susceptibility analysis.

v1: [`packages/treetime/src/timetree/inference/runner.rs`](../../packages/treetime/src/timetree/inference/runner.rs).

- `fn run_timetree` [packages/treetime/src/timetree/inference/runner.rs#L30-L79](../../packages/treetime/src/timetree/inference/runner.rs#L30-L79): computes branch distributions, then performs backward and forward propagation under the optional coalescent prior that the caller passes in.

---

## Estimation Pipeline

End-to-end timetree estimation from input to output.

v1: [`packages/treetime/src/timetree/pipeline.rs`](../../packages/treetime/src/timetree/pipeline.rs), called by the CLI command in [`packages/app-commands/src/commands/timetree/run.rs`](../../packages/app-commands/src/commands/timetree/run.rs).

- `fn run_timetree_estimation` [packages/app-commands/src/commands/timetree/run.rs#L43](../../packages/app-commands/src/commands/timetree/run.rs#L43): loads the input data, builds `TimetreeParams` from the CLI arguments, runs the pipeline and writes the outputs
- `fn run` [packages/treetime/src/timetree/pipeline.rs#L55-L141](../../packages/treetime/src/timetree/pipeline.rs#L55-L141): estimates the initial clock model, initializes the branch model, runs the pre-loop steps (ML branch-length optimization, rerooting, clock filter), runs the init phase, iterates refinement rounds with convergence monitoring, computes the final times and confidence data, and gathers the results
- `fn initialize` [packages/treetime/src/timetree/init.rs#L35](../../packages/treetime/src/timetree/init.rs#L35): seeds the expected lineage counts from the Gaussian node-time marginals at the regression rate, fits the coalescent time scale to them, runs one time inference under that prior (without a prior when the configuration requests none) and rebuilds the expected lineage counts from its posteriors

---

## Reroot

Clock-based rerooting with partition state update. Searches for the root position that maximizes temporal signal (positive clock rate, minimal regression residuals).

v1: [`packages/treetime/src/timetree/optimization/reroot.rs`](../../packages/treetime/src/timetree/optimization/reroot.rs).

- `reroot_tree()` (`#reroot_tree`) [packages/treetime/src/timetree/optimization/reroot.rs#L22-L70](../../packages/treetime/src/timetree/optimization/reroot.rs#L22-L70): performs clock regression with rerooting, applies topology changes to partitions (edge split, edge merge, inverted edges), recomputes marginal messages

---

## References

- <a id="ref-1"></a>Pearl, Judea. 1988. _Probabilistic Reasoning in Intelligent Systems: Networks of Plausible Inference._ Morgan Kaufmann. ISBN 978-0-934613-73-2. https://doi.org/10.1016/C2009-0-27609-4 [↩](#cite-1)
- <a id="ref-2"></a>Kingman, J. F. C. 1982. "The Coalescent." _Stochastic Processes and their Applications_ 13(3):235-248. https://doi.org/10.1016/0304-4149(82)90011-4 [↩](#cite-2)
- <a id="ref-3"></a>Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4(1):vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-3)
- <a id="ref-4"></a>Brent, Richard P. 1973. _Algorithms for Minimization Without Derivatives._ Prentice-Hall. ISBN 0-13-022335-2. https://maths-people.anu.edu.au/brent/pub/pub011.html [↩](#cite-4)
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

| File                                                                                                 | Algorithms                                                            |
| ---------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| [`packages/treetime/src/timetree/inference/`](../../packages/treetime/src/timetree/inference/)       | Belief propagation, branch distributions, timetree runner             |
| [`packages/treetime/src/coalescent/`](../../packages/treetime/src/coalescent/)                       | Kingman coalescent, skyline, Tc optimization                          |
| [`packages/treetime/src/timetree/optimization/`](../../packages/treetime/src/timetree/optimization/) | Polytomy, relaxed clock, reroot, clock filter                         |
| [`packages/treetime/src/timetree/convergence/`](../../packages/treetime/src/timetree/convergence/)   | Convergence monitoring, likelihood tracking, sequence change counting |
| [`packages/app-output/src/`](../../packages/app-output/src/)                                         | Output writers: tree outputs, confidence intervals, coalescent, RTT   |
| [`packages/treetime/src/timetree/round.rs`](../../packages/treetime/src/timetree/round.rs)           | EM-like iterative refinement                                          |
| [`packages/treetime/src/timetree/pipeline.rs`](../../packages/treetime/src/timetree/pipeline.rs)     | End-to-end estimation pipeline                                        |
