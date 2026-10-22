# Timetree Inference

## Initialization

- [x] Load tree from `--tree`
- [x] Date table from `--metadata`, required
- [x] Optional alignment: `--branch-length-mode` `auto` (default) keeps the input branch lengths without an alignment
- [x] Several alignment files share the branch lengths, with one relative rate per file estimated at the end of each ML branch-length pass and fixed during the dating rounds (`--partition-rates`, only with two or more files and re-estimated branch lengths, [kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md)); the clock rate of a file is the clock rate times its rate
- [x] Tips without a sequence in a file are filled as unknown; the run stops when a file misses more than a third of the tips
- [ ] Gap handling matches v0: v0 timetree with the default `--gtr infer` uses the 5-state `nuc` alphabet, where the gap is a state ([packages/legacy/treetime/treetime/wrappers.py#L56](../../packages/legacy/treetime/treetime/wrappers.py#L56)); the v1 nucleotide alphabet puts the gap in the undetermined set together with the unknown character, so a gap carries no state information ([packages/treetime/src/alphabet/alphabet.rs#L131](../../packages/treetime/src/alphabet/alphabet.rs#L131)). Named models `JC69`, `F81` and `TN93` also use the 5-state alphabet, `K80`, `HKY85` and `T92` the 4-state one ([kb/issues/M-timetree-gap-state-v0-divergence.md](../issues/M-timetree-gap-state-v0-divergence.md))
- [x] Initial node divergences always initialized
- [x] Date assignment fails when no valid dates or fewer than three leaves have valid dates
- [x] ML branch-length optimization pre-step (v0: `optimize_tree(max_iter=1)` before and after rerooting)
- [ ] Tree inference from alignment: v0 calls IQ-TREE, FastTree, or RAxML as subprocesses; v1 requires `--tree` ([kb/issues/M-timetree-tree-inference-unimplemented.md](../issues/M-timetree-tree-inference-unimplemented.md))

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
- [x] Build branch distributions from input lengths when partitions absent: Poisson for short branches, saturation-corrected Gaussian for long ones (v0 `BranchLenInterpolator` input mode)
- [x] Optional coalescent contribution calculation

## Branch-Length Modes

- [x] `auto` (default, v0 `fn _set_branch_length_mode()`): `input` without an alignment or when the longest input branch is above 0.1 substitutions per site, otherwise `marginal` (v0 chooses `joint`, which v1 does not implement); `marginal` when no input branch has a length. Resolved after the tree is read ([packages/app-commands/src/commands/timetree/branch_length_mode.rs](../../packages/app-commands/src/commands/timetree/branch_length_mode.rs))
- [x] `input` (use tree branch lengths, no sequence partitions)
- [x] `marginal` (optimize via marginal reconstruction)
  - [x] Dense marginal partition branch
  - [x] Sparse marginal partition branch (compresses sequences first)
  - [x] GTR model selection from `--model` flag (named models and inference)
- [ ] `joint` (v0 supports joint optimization)

## Time Marginal Modes

v1 has no `--time-marginal`: every round runs marginal time inference ([kb/decisions/timetree-marginal-only-time-inference.md](../decisions/timetree-marginal-only-time-inference.md)), and `--confidence` is the only switch for date intervals.

- [/] v0 `never` (default, no intervals): the v1 default; v0 uses joint most-likely times, which v1 does not implement
- [x] v0 `only-final` and `confidence-only`: `--confidence`
  - [x] Final timetree pass after loop
  - [x] Final marginal update when partitions exist
  - [x] Confidence interval extraction from node time distributions
  - [x] Lower and upper bounds in the node dates table (`dates`)
- [ ] v0 `always` (marginal times in every round, intervals without a final round): no flag, because every v1 round is marginal and `--confidence` writes the intervals

## Coalescent Models

- [x] Constant Tc (fixed from CLI `--coalescent <TC>`: substitutions per site as in v0, converted to years with the clock rate of each round, or a time such as `10y`)
- [x] Optimized Tc (`--coalescent opt`, analytic one-segment skyline solve)
  - [x] Re-estimates constant $T_c$ at the start of every loop round ([packages/treetime/src/timetree/refinement_loop.rs#L54-L61](../../packages/treetime/src/timetree/refinement_loop.rs#L54-L61))
  - [x] Pre-loop constant Tc optimization
- [/] Skyline (`--coalescent skyline`, convex Newton solve over piecewise-constant $T_c(t)$; extrapolation and quadrature contracts remain open: [kb/issues/N-coalescent-skyline-extrapolation-policy-undecided.md](../issues/N-coalescent-skyline-extrapolation-policy-undecided.md), [kb/issues/N-coalescent-skyline-quadrature-contract-undecided.md](../issues/N-coalescent-skyline-quadrature-contract-undecided.md))
  - [ ] Skyline fitted in the last round only: v0 uses a constant $T_c$ until the last round; v1 re-estimates the skyline at the start of every loop round ([packages/treetime/src/timetree/refinement_loop.rs#L54-L61](../../packages/treetime/src/timetree/refinement_loop.rs#L54-L61)) and runs no skyline step after the loop ([kb/issues/M-timetree-skyline-timing-v0-divergence.md](../issues/M-timetree-skyline-timing-v0-divergence.md))
- [x] Final marginal round after the loop for `--confidence` (v0 `only-final` and `confidence-only`): v0 runs no other post-loop timetree pass ([packages/legacy/treetime/treetime/treetime.py#L390-L395](../../packages/legacy/treetime/treetime/treetime.py#L390-L395), [packages/legacy/treetime/treetime/wrappers.py#L453-L470](../../packages/legacy/treetime/treetime/wrappers.py#L453-L470), [packages/legacy/treetime/treetime/wrappers.py#L492](../../packages/legacy/treetime/treetime/wrappers.py#L492)); v1 `fn final_marginal_round` ([packages/treetime/src/timetree/round.rs#L88-L99](../../packages/treetime/src/timetree/round.rs#L88-L99)), called from `fn refine_final_times` ([packages/treetime/src/timetree/pipeline.rs#L424-L431](../../packages/treetime/src/timetree/pipeline.rs#L424-L431))
- [x] Coalescent leaf, internal-node, and root contributions
- [x] Child-first contribution ordering matching v0
- [x] Merger rate lambda(t) = k(k-1)/(2\*Tc)
- [x] Branch counting k(t) from node times
- [ ] `--n-branches-posterior` (no flag until point-time lineage counts exist - [kb/issues/N-timetree-n-branches-posterior-unimplemented.md](../issues/N-timetree-n-branches-posterior-unimplemented.md))
- [ ] Empirical skyline (v0 `skyline_empirical()` - sliding window without optimization)
- [x] Skyline confidence bands from the full inverse Hessian ([kb/decisions/coalescent-skyline-hessian-confidence-bands.md](../decisions/coalescent-skyline-hessian-confidence-bands.md))
- [ ] Skyline plot output ([kb/issues/N-timetree-plot-unimplemented.md](../issues/N-timetree-plot-unimplemented.md))

## Polytomy Resolution

- [x] Find polytomy nodes with more than two children
- [x] Stochastic coalescent-with-mutations resolution (v0's `--stochastic-resolve`, always on in v1; [kb/decisions/timetree-stochastic-polytomy-resolution.md](../decisions/timetree-stochastic-polytomy-resolution.md))
- [x] Exact per-branch substitution counts from `edge_subs`, falling back to `round(mutation_length * L)`
- [x] Per-branch coalescent merger rate from the round's coalescent model; when no coalescent prior is requested the model is built from a constant $T_c$ estimated from the tree, replacing v0's per-polytomy window-calibrated dummy rate ([kb/decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md))
- [x] Exact hazard integration across lineage arrivals and piecewise-constant merger-rate breakpoints
- [x] Finite input, rate, event-plan, and parent-bound validation before graph mutation
- [x] `--seed` for reproducible resolution; a generated seed is logged when none is given
- [x] Remove undated single-child nodes after resolution; a dated single-child node stays with its date constraint (v0 removes it: [kb/v0-errata/timetree-polytomy-cleanup-drops-dated-single-child-node.md](../v0-errata/timetree-polytomy-cleanup-drops-dated-single-child-node.md))
- [x] Reconcile partition topology after tree change
- [x] Resolution by default as in v0; `--keep-polytomies` turns it off
- [ ] Greedy pairwise merging (v0's `--greedy-resolve`; v1 does not provide this method, and v0 deprecates it as unsuitable for large polytomies)

## Relaxed Clock

- [x] Postorder penalty coefficients (k1, k2)
- [x] Preorder optimal gamma (rate multiplier per branch)
- [x] Coupling parameter (parent-offspring rate correlation)
- [x] Slack parameter (rate deviation penalty)
- [x] Skip relaxed clock when total partition length is zero
- [x] Fit each rate against the duration between the inferred node dates times the clock rate (v0 `clock_length`)
- [x] Keep the fitted rate multiplier of every edge that polytomy resolution keeps; a new edge starts at 1.0, as in v0
- [x] Use first two values as slack and coupling, with defaults when omitted
- [x] `--relax` argument parsing (`num_args = 2` for slack and coupling)
- [ ] Substitution rates output (`substitution_rates.tsv` in v0)

## Refinement Iteration Loop

- [x] Loop stops on convergence or `--max-iter`
- [x] Convergence check uses max node-time change below `NODE_TIME_TOLERANCE_YEARS` and `n_resolved == 0`, falling back to `n_diff` when no times are comparable ([kb/decisions/timetree-convergence-on-node-times.md](../decisions/timetree-convergence-on-node-times.md))
- [x] Clock-constrained branch lengths committed after each pass and propagated by the next marginal reconstruction, damped by 0.5 ([kb/decisions/timetree-clock-constrained-profile-propagation.md](../decisions/timetree-clock-constrained-profile-propagation.md))
- [x] Coalescent lineage counts frozen before the loop for the prior; $T_c$ and the reported likelihood estimated against the live tree
- [ ] Skyline mode suppresses early convergence exit: v0 never stops early in skyline mode ([packages/legacy/treetime/treetime/treetime.py#L374](../../packages/legacy/treetime/treetime/treetime.py#L374)); v1 keeps the early exit ([kb/issues/M-timetree-skyline-timing-v0-divergence.md](../issues/M-timetree-skyline-timing-v0-divergence.md))
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
- [x] Node-time log marginal likelihood (`log_lh_time`): the backward pass's integral over node times, like v0's marginal-mode `positional_LH` (the integral of the root's marginal distribution); v0's joint-mode `positional_LH` sums branch-length likelihoods at the inferred lengths and has no v1 equivalent
- [x] Tracelog CSV output (`--output-tracelog`)
- [x] Coalescent log-likelihood in metrics (`log_lh_coal`)
- [x] Typed aggregate likelihood components (`LogLh`) with unchanged numeric CSV serialization

## Confidence Intervals

- [x] 90% maximum posterior region from marginal posteriors (`max_posterior_region()`)
- [x] `combine_in_quadrature()` quadrature sum
- [x] Inverse CDF of the posterior grid for regions with the peak on a boundary
- [x] Rate susceptibility analysis (`compute_rate_susceptibility()`)
- [x] `date_uncertainty_due_to_rate()` via probit function (erfinv)
- [x] Rate susceptibility activated via `--confidence`, with `--clock-std-dev` or the Gaussian rate standard deviation of an estimated rate (v0: `--clock-std-dev` or `--covariation`; see [kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md))
- [x] Clock rate stepped toward the maximum of the marginal likelihood of the node times in each refinement round when a coalescent prior is active - new in v1, see [timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md)
- [x] Reroot method and tip/MRCA selection passed into clock rerooting

## Output

- [x] Timetree Newick
- [x] Timetree Nexus
- [x] Tree-format topology ordering
- [x] Clock model JSON with `timetree.*` basename
- [x] Model document (`--output-model`) with one entry per alignment file, including the clock rate of each file
- [x] Per-file sequence outputs of several alignment files (FASTA, Auspice JSON, UShER MAT, Nexus, Newick with `--newick-annotations beast` or `nhx`, ancestral augur node data): `<command>.<n>.<ext>` under `--output-all`, or a `{alignment}` placeholder in the per-file flag; a per-file flag without the placeholder stops with an error
- [x] Clock regression table (`--output-rtt`, opt-in; one row per sample with the date, date source (`input` or `inferred`), root-to-tip distance and clock-filter flag the final clock model was fitted on; new in v1, v0 has no table; app runs request it for the root-to-tip plot)
- [x] Coalescent time scale TSV/CSV/JSON (`--output-coalescent-{tsv,csv,json}`; per-segment $T_c$, $N_e$, and confidence band for skyline, optimized-constant, and fixed $T_c$; TSV in the default `--output-all` set, CSV and JSON opt-in; the `skyline.pdf` plot remains, see [kb/issues/N-timetree-plot-unimplemented.md](../issues/N-timetree-plot-unimplemented.md))
- [ ] Root sequence sampled from its posterior in every reconstruction, as v0 does; v1 takes the most likely state ([kb/issues/M-timetree-root-sequence-not-sampled.md](../issues/M-timetree-root-sequence-not-sampled.md))
- [x] Ancestral sequences FASTA (`--output-sequences`, default under `--output-all`; reuses the `ancestral` marginal reconstruction, internal nodes always, tips with `--include-leaves`)
- [ ] Branch mutations table ([kb/issues/N-timetree-branch-mutations-output-unimplemented.md](../issues/N-timetree-branch-mutations-output-unimplemented.md))
- [ ] Molecular clock text output or approved replacement ([kb/issues/N-timetree-molecular-clock-text-output-undecided.md](../issues/N-timetree-molecular-clock-text-output-undecided.md))
- [ ] Sequence-evolution model text output or approved replacement ([kb/issues/N-timetree-sequence-evolution-model-text-output-undecided.md](../issues/N-timetree-sequence-evolution-model-text-output-undecided.md))
- [x] Node dates TSV (`--output-dates`, `timetree.dates.tsv`, in the default set; the columns and header of v0 `dates.tsv`, with the bounds of the 90% region with `--confidence`)
- [ ] Substitution rates TSV (v0 writes `substitution_rates.tsv` when `--relax` is used)
- [/] Auspice v2 JSON (schema-valid substitutions, dates, divergence, outlier state, sequences, and genome annotations; no Augur golden master checks the document: [kb/issues/N-io-auspice-output-missing-augur-golden-master.md](../issues/N-io-auspice-output-missing-augur-golden-master.md))
- [ ] Outliers TSV (v0 writes `outliers.tsv`)
- [ ] Tracelog run (v0 `tracelog_run()` with detailed per-iteration state)
- [ ] Plotting (`--plot-tree`, `--plot-rtt`, `--tip-labels`, `--no-tip-labels` - no flags until the plots exist - [kb/issues/N-timetree-plot-unimplemented.md](../issues/N-timetree-plot-unimplemented.md))

## CLI Options (Parsed but Not Wired)

- [x] `--clock-std-dev` (rate susceptibility; only with `--confidence`)
- [x] `--confidence` (final marginal round, rate susceptibility)
- [ ] `--covariation`, `--tip-slack`: removed in v1 ([kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md))
- [x] `--reroot` / `--outgroup` / `--keep-root`
- [x] `--model` (model selection: named models and inference)
- [x] `--model-params` (`key=value` parameters of a named model, as v0 `--gtr-params`; an unknown key or a malformed item is an error)
- [x] `--zero-based` (mutation positions in the comments of Newick and Nexus trees)
- [x] `--include-leaves` / `--impute-missing-data` / `--reconstruct-tip-states` (drive the reconstructed FASTA: emit tips, impute ambiguous tip states, alias enabling both; see [kb/decisions/ancestral-marginal-tip-reconstruction-and-imputation.md](../decisions/ancestral-marginal-tip-reconstruction-and-imputation.md))
- [x] `--report-ambiguous` (keep mutations from or to `N`/`X` in written mutation lists; omitted by default, as in v0; see [kb/decisions/ancestral-marginal-tip-reconstruction-and-imputation.md](../decisions/ancestral-marginal-tip-reconstruction-and-imputation.md))
- [x] `--seed` (see [kb/decisions/cli-seed-on-commands-with-random-steps.md](../decisions/cli-seed-on-commands-with-random-steps.md))
- [x] `--gen-per-year` (generations per year, default 50.0, only with `--coalescent`; reports effective population size `N_e = Tc * gen_per_year` to the log for the constant, opt, and skyline modes, and to the coalescent output files above)
- [ ] `--aln` legacy option
