# Timetree Inference

## Initialization

- [x] Load tree from `--tree`
- [x] Optional date table from `--dates`
- [x] Optional alignment when `--branch-length-mode=marginal`
- [x] Alphabet always maps gap to unknown profile (matching v0 `nuc_nogap` behavior)
- [x] Initial node divergences always initialized
- [x] Date assignment fails when no valid dates or fewer than three leaves have valid dates
- [x] ML branch-length optimization pre-step (v0: `optimize_tree(max_iter=1)` before and after rerooting)
- [ ] Tree inference from alignment path is `todo!`

## Time Distribution Propagation

- [x] Backward pass (leaf dates to root via convolution)
- [x] Forward pass (root to leaves, refine distributions)
  - [x] Refines every node whose date is not exact, leaves included (uncertain, ranged, and missing dates)
  - [x] Leaves a node given an exact date at that date, unrefined and unclamped
  - [x] Keeps a given date the rest of the tree contradicts, rather than leaving the node undated
  - [x] Cuts the parent to the times that can reach the node before convolving
- [x] Backward pass lifts the fixed input date constraint back into the time distribution each round
- [x] Bad branch exclusion (outliers, dateless leaves)
- [x] Build branch distributions from partitions when present
- [x] Build point branch distributions from input lengths when partitions absent
- [x] Optional coalescent contribution calculation

## Branch-Length Modes

- [x] `input` (use tree branch lengths, no sequence partitions, point distributions)
- [x] `marginal` (optimize via marginal reconstruction)
  - [x] Dense marginal partition branch
  - [x] Sparse marginal partition branch (compresses sequences first)
  - [x] GTR model selection from `--gtr` flag (named models and inference)
- [ ] `joint` (v0 supports joint optimization)

## Time Marginal Modes

- [x] `never` (joint most-likely times)
- [ ] `always` (marginal every round - parsed but not wired)
- [x] `only-final` (marginal last round for confidence)
  - [x] Final timetree pass after loop
  - [x] Final marginal update when partitions exist
  - [x] Confidence interval extraction from node time distributions
  - [x] Confidence TSV output
- [x] `confidence-only` (v0 variant, implemented as `--confidence` promoting `Never` to `OnlyFinal`)

## Coalescent Models

- [x] Constant Tc (fixed from CLI `--coalescent`)
- [x] Optimized Tc (`--coalescent-opt`, analytic one-segment skyline solve)
  - [x] Re-optimizes constant Tc inside loop for iterations i >= 2
  - [x] Pre-loop constant Tc optimization
- [/] Skyline (`--coalescent-skyline`, convex Newton solve over piecewise-constant $T_c(t)$; extrapolation, quadrature, and internal boundary-validation contracts remain open: [kb/issues/N-coalescent-skyline-extrapolation-policy-undecided.md](../issues/N-coalescent-skyline-extrapolation-policy-undecided.md), [kb/issues/N-coalescent-skyline-quadrature-contract-undecided.md](../issues/N-coalescent-skyline-quadrature-contract-undecided.md), [kb/issues/N-coalescent-skyline-grid-validation-incomplete.md](../issues/N-coalescent-skyline-grid-validation-incomplete.md))
  - [x] Final skyline re-optimization after refinement loop
  - [x] Extra final timetree pass unless `--time-marginal=only-final`
- [x] Coalescent leaf, internal-node, and root contributions
- [x] Child-first contribution ordering matching v0
- [x] Merger rate lambda(t) = k(k-1)/(2\*Tc)
- [x] Branch counting k(t) from node times
- [ ] `--n-branches-posterior` (parsed, returns explicit error - [kb/issues/N-timetree-n-branches-posterior-unimplemented.md](../issues/N-timetree-n-branches-posterior-unimplemented.md))
- [ ] Empirical skyline (v0 `skyline_empirical()` - sliding window without optimization)
- [ ] Skyline confidence intervals (v0 computes via second derivatives)
- [ ] Skyline plot output ([kb/issues/N-timetree-missing-skyline-output.md](../issues/N-timetree-missing-skyline-output.md))

## Polytomy Resolution

- [x] Find polytomy nodes with more than two children
- [x] Stochastic coalescent-with-mutations resolution (v0's `--stochastic-resolve`, always on in v1; [kb/decisions/timetree-stochastic-polytomy-resolution.md](../decisions/timetree-stochastic-polytomy-resolution.md))
- [x] Exact per-branch substitution counts from `edge_subs`, falling back to `round(mutation_length * L)`
- [x] Per-branch coalescent merger rate from the round's coalescent model; when no coalescent prior is requested the model is built from a constant $T_c$ estimated from the tree, replacing v0's per-polytomy window-calibrated dummy rate ([kb/decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md))
- [x] Exact hazard integration across lineage arrivals and piecewise-constant merger-rate breakpoints
- [x] Finite input, rate, event-plan, and parent-bound validation before graph mutation
- [x] `--seed` for reproducible resolution; a generated seed is logged when none is given
- [x] Remove obsolete single-child nodes after resolution
- [x] Reconcile partition topology after tree change
- [/] `--keep-polytomies` (parsed but never read - [kb/issues/N-timetree-unused-cli-flags.md](../issues/N-timetree-unused-cli-flags.md))
- [ ] Greedy pairwise merging (v0's `--greedy-resolve`; v1 does not provide this method, and v0 deprecates it as unsuitable for large polytomies)

## Relaxed Clock

- [x] Postorder penalty coefficients (k1, k2)
- [x] Preorder optimal gamma (rate multiplier per branch)
- [x] Coupling parameter (parent-offspring rate correlation)
- [x] Slack parameter (rate deviation penalty)
- [x] Skip relaxed clock when total partition length is zero
- [x] Use first two values as slack and coupling, with defaults when omitted
- [x] `--relax` argument parsing (`num_args = 2` for slack and coupling)
- [ ] Substitution rates output (`substitution_rates.tsv` in v0)

## Refinement Iteration Loop

- [x] Loop stops on convergence or `--max-iter`
- [x] Convergence check uses max node-time change below `NODE_TIME_TOLERANCE_YEARS` and `n_resolved == 0`, falling back to `n_diff` when no times are comparable ([kb/decisions/timetree-convergence-on-node-times.md](../decisions/timetree-convergence-on-node-times.md))
- [x] Clock-constrained branch lengths committed after each pass and propagated by the next marginal reconstruction, damped by 0.5 ([kb/decisions/timetree-clock-constrained-profile-propagation.md](../decisions/timetree-clock-constrained-profile-propagation.md))
- [x] Coalescent lineage counts frozen before the loop for the prior; $T_c$ and the reported likelihood estimated against the live tree
- [x] Skyline mode suppresses early convergence exit
- [x] Relaxed clock application
- [x] Polytomy resolution
- [x] Dirty-tree-aware reconstruction ordering
  - [x] If topology changed, reconcile partitions, propagate internal bad-branch state, rerun marginal reconstruction, establish node times without coalescent, then rerun timetree with coalescent
  - [x] If topology did not change, rerun marginal first and timetree second
- [x] Re-estimate clock model after every iteration without rerooting
- [x] Optional tracelog write on each recorded iteration

## Convergence Tracking

- [x] Node-time movement per round (`max_time_change`, `rms_time_change`) - new in v1
- [x] Sequence change count (n_diff)
- [x] Polytomies resolved count (n_resolved)
- [x] Sequence log-likelihood (`log_lh_seq`)
- [x] Positional log-likelihood (`log_lh_pos`)
- [x] Tracelog CSV output (`--tracelog`)
- [x] Coalescent log-likelihood in metrics (`log_lh_coal`)
- [x] Typed aggregate likelihood components (`LogLh`) with unchanged numeric CSV serialization

## Confidence Intervals

- [x] 90% HPD from marginal posteriors
- [x] `combine_confidence()` quadrature sum
- [x] `Distribution::quantile()` inverse CDF
- [x] Rate susceptibility analysis (`compute_rate_susceptibility()`)
- [x] `date_uncertainty_due_to_rate()` via probit function (erfinv)
- [x] Rate susceptibility activated via `--confidence` with `--covariation` or `--clock-std-dev`
- [x] `--confidence` promotion of `time_marginal` from `never` to `only-final`
- [x] `--covariation` wired into timetree clock regression
- [x] Reroot method and tip/MRCA selection passed into clock rerooting

## Output

- [x] Timetree Newick
- [x] Timetree Nexus
- [x] Tree-format topology ordering
- [x] Clock model JSON with `timetree.*` basename
- [x] Confidence TSV
- [ ] Ancestral sequences FASTA ([kb/issues/N-timetree-ancestral-sequences-output-unimplemented.md](../issues/N-timetree-ancestral-sequences-output-unimplemented.md))
- [ ] Branch mutations table ([kb/issues/N-timetree-branch-mutations-output-unimplemented.md](../issues/N-timetree-branch-mutations-output-unimplemented.md))
- [ ] Molecular clock text output or approved replacement ([kb/issues/N-timetree-molecular-clock-text-output-undecided.md](../issues/N-timetree-molecular-clock-text-output-undecided.md))
- [ ] Sequence-evolution model text output or approved replacement ([kb/issues/N-timetree-sequence-evolution-model-text-output-undecided.md](../issues/N-timetree-sequence-evolution-model-text-output-undecided.md))
- [ ] Node dates TSV (`write_node_dates()` is `todo!()` - [kb/issues/N-timetree-node-dates-output-unimplemented.md](../issues/N-timetree-node-dates-output-unimplemented.md))
- [ ] Substitution rates TSV (v0 writes `substitution_rates.tsv` when `--relax` is used)
- [/] Auspice v2 JSON (schema-valid substitutions, dates, divergence, outlier state, sequences, and genome annotations; entropy perturbs the Shannon definition and inference metadata is incomplete: [kb/issues/M-io-auspice-entropy-perturbs-shannon-definition.md](../issues/M-io-auspice-entropy-perturbs-shannon-definition.md), [kb/issues/M-timetree-tree-output-inference-metadata-incomplete.md](../issues/M-timetree-tree-output-inference-metadata-incomplete.md))
- [ ] Outliers TSV (v0 writes `outliers.tsv`)
- [ ] Tracelog run (v0 `tracelog_run()` with detailed per-iteration state)
- [ ] Plotting (`--plot-tree`, `--plot-rtt` - parsed, return explicit error - [kb/issues/N-timetree-plot-unimplemented.md](../issues/N-timetree-plot-unimplemented.md))

## CLI Options (Parsed but Not Wired)

- [x] `--clock-std-dev` (rate susceptibility)
- [x] `--confidence` (time_marginal promotion, rate susceptibility)
- [x] `--covariation` (GLS clock regression, rate susceptibility)
- [x] `--tip-slack` (covariation variance computation)
- [x] `--reroot` / `--reroot-tips`
- [ ] `--n-iqd`
- [ ] `--tip-labels` / `--no-tip-labels`
- [x] `--gtr` (model selection: named models and inference)
- [ ] `--gtr-params` (parsed but not wired)
- [ ] `--method-anc`
- [ ] `--keep-overhangs` (gap handling not implemented)
- [ ] `--zero-based` ([kb/issues/M-core-mutation-representation-and-format-projection-inconsistent.md](../issues/M-core-mutation-representation-and-format-projection-inconsistent.md))
- [ ] `--reconstruct-tip-states`
- [ ] `--report-ambiguous`
- [ ] `--seed`
- [x] `--gen-per-year` (generations per year, default 50.0; reports effective population size `N_e = Tc * gen_per_year` to the log for the constant, opt, and skyline modes; file output tracked with skyline plot output above)
- [ ] `--aln` legacy option
