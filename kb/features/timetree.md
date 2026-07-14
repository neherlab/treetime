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
- [x] Bad branch exclusion (outliers, dateless leaves)
- [x] Build branch distributions from partitions when present
- [x] Build point branch distributions from input lengths when partitions absent
- [/] Optional coalescent contribution calculation ([kb/issues/M-timetree-coalescent-missing-leaf-and-root-contributions.md](../issues/M-timetree-coalescent-missing-leaf-and-root-contributions.md))

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
- [/] Optimized Tc (`--coalescent-opt`, Brent's method, log-space bracket [-20, 2]; the command can report success without retaining an optimized parameter: [kb/issues/M-timetree-constant-tc-success-without-parameter.md](../issues/M-timetree-constant-tc-success-without-parameter.md))
  - [x] Re-optimizes constant Tc inside loop for iterations i >= 2
  - [/] No special pre-loop setup unless skyline also active
- [/] Skyline (`--coalescent-skyline`, Nelder-Mead, piecewise linear $T_c(t)$; simplex, extrapolation, quadrature, and grid-validation contracts remain open: [kb/issues/N-coalescent-skyline-simplex-initialization-undecided.md](../issues/N-coalescent-skyline-simplex-initialization-undecided.md), [kb/issues/N-coalescent-skyline-extrapolation-policy-undecided.md](../issues/N-coalescent-skyline-extrapolation-policy-undecided.md), [kb/issues/N-coalescent-skyline-quadrature-contract-undecided.md](../issues/N-coalescent-skyline-quadrature-contract-undecided.md), [kb/issues/N-coalescent-skyline-grid-validation-incomplete.md](../issues/N-coalescent-skyline-grid-validation-incomplete.md))
  - [x] Pre-loop constant Tc optimization
  - [x] Warning fallback to Tc = 1.0 on failure or non-convergence
  - [x] Final skyline re-optimization after refinement loop
  - [x] Extra final timetree pass unless `--time-marginal=only-final`
- [/] Coalescent contribution per node (leaf/root terms are incomplete - [kb/issues/M-timetree-coalescent-missing-leaf-and-root-contributions.md](../issues/M-timetree-coalescent-missing-leaf-and-root-contributions.md))
- [/] Different multiplication ordering ([kb/decisions/coalescent-multiplication-ordering.md](../decisions/coalescent-multiplication-ordering.md), [kb/issues/M-timetree-coalescent-multiplication-ordering-diverges-from-v0.md](../issues/M-timetree-coalescent-multiplication-ordering-diverges-from-v0.md))
- [x] Merger rate lambda(t) = k(k-1)/(2\*Tc)
- [x] Branch counting k(t) from node times
- [ ] `--n-branches-posterior` (parsed, returns explicit error - [kb/issues/N-timetree-n-branches-posterior-unimplemented.md](../issues/N-timetree-n-branches-posterior-unimplemented.md))
- [ ] Empirical skyline (v0 `skyline_empirical()` - sliding window without optimization)
- [ ] Skyline confidence intervals (v0 computes via second derivatives)
- [ ] Skyline plot output ([kb/issues/N-timetree-missing-skyline-output.md](../issues/N-timetree-missing-skyline-output.md))

## Polytomy Resolution

- [x] Find polytomy nodes with more than two children
- [x] Greedy pairwise merging (likelihood gain threshold 0.05)
- [x] Per-pair Brent optimization for merge time
- [x] Remove obsolete single-child nodes after resolution
- [x] Reconcile partition topology after tree change
- [/] `--keep-polytomies` (parsed but never read - [kb/issues/N-timetree-dead-cli-flags.md](../issues/N-timetree-dead-cli-flags.md))
- [ ] Stochastic resolution (`--stochastic-resolve` in v0 - [kb/issues/N-timetree-stochastic-polytomy-unimplemented.md](../issues/N-timetree-stochastic-polytomy-unimplemented.md))

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
- [x] Convergence check uses `n_diff == 0 && n_resolved == 0`
- [x] Skyline mode suppresses early convergence exit
- [x] Relaxed clock application
- [x] Polytomy resolution
- [/] Dirty-tree-aware reconstruction ordering (the sequence of passes is implemented, but coalescent events can remain incomplete after topology changes: [kb/issues/H-timetree-coalescent-events-incomplete-after-topology-change.md](../issues/H-timetree-coalescent-events-incomplete-after-topology-change.md))
  - [x] If topology changed, rerun timetree first and marginal second
  - [x] If topology did not change, rerun marginal first and timetree second
- [x] Re-estimate clock model after every iteration without rerooting
- [x] Optional tracelog write on each recorded iteration

## Convergence Tracking

- [x] Sequence change count (n_diff)
- [x] Polytomies resolved count (n_resolved)
- [x] Sequence likelihood (lh_seq)
- [x] Positional likelihood (lh_pos)
- [x] Tracelog CSV output (`--tracelog`)
- [x] Coalescent likelihood in metrics (lh_coal)

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
- [/] Tree-format topology ordering ([kb/issues/M-io-tree-backed-output-order-inconsistent.md](../issues/M-io-tree-backed-output-order-inconsistent.md))
- [x] Clock model JSON with `timetree.*` basename
- [x] Confidence TSV
- [ ] Node dates TSV (`write_node_dates()` is `todo!()` - [kb/issues/N-timetree-node-dates-output-unimplemented.md](../issues/N-timetree-node-dates-output-unimplemented.md))
- [ ] Substitution rates TSV (v0 writes `substitution_rates.tsv` when `--relax` is used)
- [/] Auspice JSON (v1 writes substitutions, `num_date`, `div`, and `bad_branch`; required `meta.updated` and remaining inference metadata are incomplete: [kb/issues/H-io-auspice-v2-required-updated-missing.md](../issues/H-io-auspice-v2-required-updated-missing.md), [kb/issues/M-timetree-tree-output-inference-metadata-incomplete.md](../issues/M-timetree-tree-output-inference-metadata-incomplete.md))
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
- [ ] `--gen-per-year` (generations per year for N_e estimation)
- [ ] `--aln` legacy option
