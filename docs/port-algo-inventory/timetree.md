# Timetree Inference Algorithms

[Back to index](_index.md)

## Belief Propagation

| Property    | Value                                                                                                                                                                                                                                                                                            |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Type        | Well-known                                                                                                                                                                                                                                                                                       |
| v1 Location | [`packages/treetime/src/commands/timetree/inference/backward_pass.rs`](../../packages/treetime/src/commands/timetree/inference/backward_pass.rs), [`packages/treetime/src/commands/timetree/inference/forward_pass.rs`](../../packages/treetime/src/commands/timetree/inference/forward_pass.rs) |
| v0 Location | [`packages/legacy/treetime/treetime/node_interpolator.py`](../../packages/legacy/treetime/treetime/node_interpolator.py)                                                                                                                                                                         |
| Reference   | Pearl, J. (1988). "Probabilistic Reasoning in Intelligent Systems." Morgan Kaufmann                                                                                                                                                                                                              |

Two-pass message passing for time inference: backward (convolve + multiply child distributions), forward (divide + convolve parent distributions).

- `propagate_distributions_backward()` (`#propagate_distributions_backward`) [packages/treetime/src/commands/timetree/inference/backward_pass.rs#L17-L31](../../packages/treetime/src/commands/timetree/inference/backward_pass.rs#L17-L31): skips bad branches (outlier and dateless leaves) so they do not constrain parent time
- `propagate_distributions_forward()` (`#propagate_distributions_forward`) [packages/treetime/src/commands/timetree/inference/forward_pass.rs#L11-L22](../../packages/treetime/src/commands/timetree/inference/forward_pass.rs#L11-L22): preserves internal node times when forward pass yields `None`

---

## Branch Length Distributions

| Property    | Value                                                                                                                                                                  |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Custom (TreeTime-specific)                                                                                                                                             |
| v1 Location | [`packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs`](../../packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs) |

Computes per-edge time distributions from partition contributions, clock rate, and gamma rate multiplier.

- `compute_branch_length_distribution()` (`#compute_branch_length_distribution`) [packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L31-L63](../../packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L31-L63): converts branch length grid to time grid using `effective_clock_rate = clock_rate * gamma`, where `gamma > 1` means faster evolution (shorter time for same substitutions)

---

## Kingman Coalescent

| Property    | Value                                                                                                                        |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                   |
| v1 Location | [`packages/treetime/src/commands/timetree/coalescent/`](../../packages/treetime/src/commands/timetree/coalescent/) (9 files) |
| v0 Location | [`packages/legacy/treetime/treetime/merger_models.py`](../../packages/legacy/treetime/treetime/merger_models.py)             |
| Reference   | Kingman, J.F.C. (1982). "The coalescent." Stochastic Processes and Applications, 13(3):235-248                               |

Computes lineage dynamics k(t), integral merger rate I(t), per-node survival/merger contributions.

- `compute_coalescent_contributions()` (`#compute_coalescent_contributions`) [packages/treetime/src/commands/timetree/coalescent/coalescent.rs#L61-L79](../../packages/treetime/src/commands/timetree/coalescent/coalescent.rs#L61-L79): orchestrator that chains event collection, lineage counting, integral merger rate, and node contributions
- `compute_node_contributions()` (`#compute_node_contributions`) [packages/treetime/src/commands/timetree/coalescent/contributions.rs#L26-L71](../../packages/treetime/src/commands/timetree/coalescent/contributions.rs#L26-L71): leaf contribution = -I(t) (survival), internal = multiplicity \* (I(t) - log(lambda(t))) (merger density)

**v1 initialization strategy**: first timetree pass runs without coalescent to establish node time distributions via backward+forward pass. Coalescent contributions are computed from these established times on the second pass.

---

## Tc Optimization

| Property    | Value                                                                                                                                                  |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Type        | Well-known optimization                                                                                                                                |
| v1 Location | [`packages/treetime/src/commands/timetree/coalescent/optimize_tc.rs#L44`](../../packages/treetime/src/commands/timetree/coalescent/optimize_tc.rs#L44) |
| Reference   | Brent 1973                                                                                                                                             |

Optimizes coalescent time scale in log space [-20, 2] using Brent's method.

- `optimize_tc()` (`#optimize_tc`) [packages/treetime/src/commands/timetree/coalescent/optimize_tc.rs#L44-L84](../../packages/treetime/src/commands/timetree/coalescent/optimize_tc.rs#L44-L84): precomputes lineage counts and node branch data, then minimizes negative total coalescent likelihood
- Nodes with undetermined branch length are skipped with a warning

**v1 loop integration**: Tc is re-optimized each iteration (from iteration 2 onward) using constant Tc. In skyline mode, constant Tc is used during loop iterations; full skyline fit is deferred to post-convergence.

---

## Skyline Coalescent

| Property    | Value                                                                                                                                          |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known (Nelder-Mead) + domain-specific (skyline)                                                                                           |
| v1 Location | [`packages/treetime/src/commands/timetree/coalescent/skyline.rs#L68`](../../packages/treetime/src/commands/timetree/coalescent/skyline.rs#L68) |
| v0 Location | [`packages/legacy/treetime/treetime/merger_models.py#L281`](../../packages/legacy/treetime/treetime/merger_models.py#L281)                     |
| Reference   | Strimmer & Pybus (2001). "Exploring the demographic history." Mol Biol Evol, 18(12):2298-2305                                                  |
| Reference   | Minin et al. (2008). "Smooth skyride through a rough skyline." Mol Biol Evol, 25(7):1459-1471                                                  |

Piecewise-varying Tc estimation with smoothness regularization.

- `optimize_skyline()` (`#optimize_skyline`) [packages/treetime/src/commands/timetree/coalescent/skyline.rs#L68-L157](../../packages/treetime/src/commands/timetree/coalescent/skyline.rs#L68-L157): minimizes negative log-likelihood + stiffness penalty + boundary penalty using Nelder-Mead
- Skyline is re-optimized after iterations converge with stabilized node times

**Modification**: v1 uses Nelder-Mead; v0 uses SLSQP.

---

## Relaxed Clock

| Property    | Value                                                                                                                                                          |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                                                     |
| v1 Location | [`packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs#L25`](../../packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs#L25) |
| Reference   | Thorne et al. (1998). "Estimating the rate of evolution of the rate of molecular evolution." Mol Biol Evol, 15(12):1647-1657                                   |

Two-pass quadratic penalty optimization for rate variation.

- `apply_relaxed_clock()` (`#apply_relaxed_clock`) [packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs#L25-L125](../../packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs#L25-L125): postorder pass computes k1/k2 penalty coefficients, preorder pass computes optimal gamma (rate multiplier) per branch. Division by zero guarded.
- `one_mutation` computed as sum of sequence lengths across all partitions

---

## Polytomy Resolution

| Property    | Value                                                                                                                                                |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Custom (TreeTime-specific)                                                                                                                           |
| v1 Location | [`packages/treetime/src/commands/timetree/optimization/polytomy.rs#L27`](../../packages/treetime/src/commands/timetree/optimization/polytomy.rs#L27) |

Greedy pairwise merging with Brent optimization per pair.

- `resolve_polytomies()` (`#resolve_polytomies`) [packages/treetime/src/commands/timetree/optimization/polytomy.rs#L27-L32](../../packages/treetime/src/commands/timetree/optimization/polytomy.rs#L27-L32): entry point with default threshold (0.05)
- `prepare_tree_after_topology_change()` (`#prepare_tree_after_topology_change`) [packages/treetime/src/commands/timetree/optimization/polytomy.rs#L432-L454](../../packages/treetime/src/commands/timetree/optimization/polytomy.rs#L432-L454): clears cached distributions on internal nodes only; leaf date constraints and bad_branch flags are preserved
- After resolution, partition data is reconciled via `reconcile_topology()` to ensure entries exist for new nodes/edges

**Note**: Stochastic resolution (v0's `generate_subtree()` (`#generate_subtree`)) not ported. See [unimplemented](unimplemented.md#stochastic-polytomy-resolution).

---

## Clock Filter / Outlier Detection

| Property    | Value                                                                                                                                                |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known (IQD-based)                                                                                                                               |
| v1 Location | [`packages/treetime/src/commands/timetree/optimization/clock_filter.rs`](../../packages/treetime/src/commands/timetree/optimization/clock_filter.rs) |

Marks outlier leaves as bad branches and propagates upward.

- `apply_outlier_bad_branches()` (`#apply_outlier_bad_branches`) [packages/treetime/src/commands/timetree/optimization/clock_filter.rs#L75-L95](../../packages/treetime/src/commands/timetree/optimization/clock_filter.rs#L75-L95): sets `bad_branch=true` on outlier leaves, then postorder propagation marks internal nodes bad only when all children are bad

---

## Convergence Monitoring

| Property    | Value                                                                                                                          |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------ |
| Type        | Custom (TreeTime-specific)                                                                                                     |
| v1 Location | [`packages/treetime/src/commands/timetree/convergence/`](../../packages/treetime/src/commands/timetree/convergence/) (3 files) |

Tracks optimization loop convergence via sequence change counts and likelihood components.

- `TimetreeOptimizer` (`#TimetreeOptimizer`) [packages/treetime/src/commands/timetree/convergence/metrics.rs#L16-L108](../../packages/treetime/src/commands/timetree/convergence/metrics.rs#L16-L108): iteration controller that records `ConvergenceMetrics` per round and stops when converged or max iterations reached. Convergence criterion: `n_diff == 0 && n_resolved == 0`. Supports CSV tracelog output. Convergence can be suppressed (for skyline mode where constant Tc is used during iterations).
- `count_sequence_changes()` (`#count_sequence_changes`) [packages/treetime/src/commands/timetree/convergence/sequence_changes.rs#L19-L47](../../packages/treetime/src/commands/timetree/convergence/sequence_changes.rs#L19-L47): compares ancestral state snapshots position-by-position across all internal nodes and partitions
- `capture_ancestral_states()` (`#capture_ancestral_states`) [packages/treetime/src/commands/timetree/convergence/sequence_changes.rs#L50-L70](../../packages/treetime/src/commands/timetree/convergence/sequence_changes.rs#L50-L70): snapshots reconstructed sequences for all internal nodes. Captured before polytomy resolution to avoid inflating n_diff with newly created nodes.

### Likelihood Components

- `compute_sequence_likelihood()` (`#compute_sequence_likelihood`) [packages/treetime/src/commands/timetree/convergence/likelihood.rs#L11-L25](../../packages/treetime/src/commands/timetree/convergence/likelihood.rs#L11-L25): sum of per-partition root log-likelihoods from marginal reconstruction
- `compute_positional_likelihood()` (`#compute_positional_likelihood`) [packages/treetime/src/commands/timetree/convergence/likelihood.rs#L35-L76](../../packages/treetime/src/commands/timetree/convergence/likelihood.rs#L35-L76): sum of log-probabilities of branch length distributions evaluated at inferred time durations. **v1-specific metric** - v0's `positional_LH` sums node-level marginal log-likelihoods from the forward pass. Both trend in the same direction during convergence but produce different numerical values.
- `compute_coalescent_likelihood()` (`#compute_coalescent_likelihood`) [packages/treetime/src/commands/timetree/convergence/likelihood.rs#L81-L83](../../packages/treetime/src/commands/timetree/convergence/likelihood.rs#L81-L83): stub, returns `None`

---

## Iterative EM-like Refinement

| Property    | Value                                                                                                                          |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------ |
| Type        | Well-known (EM)                                                                                                                |
| v1 Location | [`packages/treetime/src/commands/timetree/refinement.rs#L21`](../../packages/treetime/src/commands/timetree/refinement.rs#L21) |
| Reference   | Sagulenko et al. (2018). "TreeTime." Virus Evolution, 4(1):vex042, Section 2.4                                                 |

Alternates sequence reconstruction (E-step) and time inference (M-step).

- `run_refinement_iteration()` (`#run_refinement_iteration`) [packages/treetime/src/commands/timetree/refinement.rs#L21-L106](../../packages/treetime/src/commands/timetree/refinement.rs#L21-L106): per-iteration logic: relaxed clock, polytomy resolution, ancestral reconstruction, timetree inference, clock re-estimation. Captures ancestral state snapshots before polytomy resolution.

---

## Confidence Intervals

| Property    | Value                                                                                                                                |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| Type        | Well-known (quantile-based)                                                                                                          |
| v1 Location | [`packages/treetime/src/commands/timetree/output/confidence.rs`](../../packages/treetime/src/commands/timetree/output/confidence.rs) |

Extracts 95% confidence intervals from marginal posterior time distributions.

- `extract_confidence_intervals()` (`#extract_confidence_intervals`) [packages/treetime/src/commands/timetree/output/confidence.rs#L46-L71](../../packages/treetime/src/commands/timetree/output/confidence.rs#L46-L71): computes 95% CI for each node using quantiles at 0.025 and 0.975. Delta distributions yield identity interval `[date, date]`.
- `combine_confidence()` (`#combine_confidence`) [packages/treetime/src/commands/timetree/output/confidence.rs#L88-L105](../../packages/treetime/src/commands/timetree/output/confidence.rs#L88-L105): combines independent CI contributions via quadrature sum
- `compute_rate_susceptibility()` (`#compute_rate_susceptibility`) [packages/treetime/src/commands/timetree/output/confidence.rs#L34-L40](../../packages/treetime/src/commands/timetree/output/confidence.rs#L34-L40): stub for clock rate uncertainty propagation

---

## Timetree Runner

| Property    | Value                                                                                                                              |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Pipeline orchestration                                                                                                             |
| v1 Location | [`packages/treetime/src/commands/timetree/inference/runner.rs`](../../packages/treetime/src/commands/timetree/inference/runner.rs) |

Single timetree inference pass: branch distributions, coalescent contributions, backward pass, forward pass.

- `run_timetree()` (`#run_timetree`) [packages/treetime/src/commands/timetree/inference/runner.rs#L23-L74](../../packages/treetime/src/commands/timetree/inference/runner.rs#L23-L74): computes branch distributions (from partitions or input lengths), optional coalescent contributions, then backward+forward propagation. Validates clock rate is positive.

---

## Estimation Pipeline

| Property    | Value                                                                                                    |
| ----------- | -------------------------------------------------------------------------------------------------------- |
| Type        | Pipeline orchestration                                                                                   |
| v1 Location | [`packages/treetime/src/commands/timetree/run.rs`](../../packages/treetime/src/commands/timetree/run.rs) |

End-to-end timetree estimation from input to output.

- `run_timetree_estimation()` (`#run_timetree_estimation`) [packages/treetime/src/commands/timetree/run.rs#L32-L253](../../packages/treetime/src/commands/timetree/run.rs#L32-L253): loads data, estimates clock model, reroots, runs clock filter, initializes partitions, runs initial timetree pass (without coalescent), optimizes Tc, runs second pass (with coalescent), iterates refinement with convergence monitoring, post-convergence skyline optimization, final marginal reconstruction for CI extraction, writes outputs.

---

## Reroot

| Property    | Value                                                                                                                                            |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| Type        | Well-known                                                                                                                                       |
| v1 Location | [`packages/treetime/src/commands/timetree/optimization/reroot.rs#L20`](../../packages/treetime/src/commands/timetree/optimization/reroot.rs#L20) |

Clock-based rerooting with partition state update.

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
