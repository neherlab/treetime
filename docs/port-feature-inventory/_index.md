# Feature Inventory - TreeTime v0/v1

> **Status**: reviewed against source code
>
> **Valid for commit**: `404d8b77` on 2026-03-10
>
> **Scope**: v0/v1 feature parity tracking. Combines domain taxonomy with CLI wiring status.

## Legend

Uses [Obsidian checkbox statuses](https://publish.obsidian.md/tasks/Getting+Started/Statuses):

- `[x]` implemented in v1 with v0 parity (or v1-only feature)
- `[/]` partial implementation, stubbed, or only partly wired
- `[ ]` missing, parsed but not wired, or documented but unimplemented

## Command Map

| Command     | Status | Notes                                      |
| ----------- | ------ | ------------------------------------------ |
| `ancestral` | [x]    | Parsimony and marginal reconstruction      |
| `clock`     | [x]    | Regression and rerooting                   |
| `timetree`  | [x]    | Full inference pipeline                    |
| `optimize`  | [x]    | v1-only branch-length optimization         |
| `prune`     | [x]    | v1-only tree pruning                       |
| `mugration` | [x]    | Marginal reconstruction with iterative GTR |
| `homoplasy` | [ ]    | Unimplemented                              |

## 1. Ancestral Reconstruction

### Fitch Parsimony (Sparse)

- [x] Backward pass (intersection/union of child state sets)
- [x] Forward pass (resolve ambiguities top-down)
- [x] Indel tracking (insertions, deletions)
- [x] Composition tracking (character counts)
- [x] Deterministic root state resolution ([intentional change](../port-intentional-changes/ancestral-fitch-deterministic-root-state.md))
- [x] Forward cleanup drops stored full sequences from non-root internal nodes

### Marginal Reconstruction (Dense)

- [x] Backward pass (log-space message multiplication)
- [x] Forward pass (parent-child message propagation)
- [x] Root equilibrium frequency weighting
- [x] Log-space normalization (logsumexp)
- [x] expQt matrix propagation
- [x] `initialize_marginal()` bootstrap with dummy JC69
- [x] Pre-inference marginal pass to populate profiles
- [x] GTR selection after profiles exist
- [x] Second marginal pass after final GTR assignment
- [ ] Profile sampling (`sample_from_prof` - v0 supports stochastic sampling, v1 argmax only)

### Marginal Reconstruction (Sparse)

- [x] Backward pass (variable positions only)
- [x] Forward pass (variable positions only)
- [x] Fixed-position per-character profiles
- [x] Mutation composition from edge substitutions
- [x] Message combine with variability threshold (EPS=1e-4)
- [x] Sparse compression before GTR selection
- [x] Per-partition GTR assignment after dummy JC69 bootstrap
- [x] Single `update_marginal()` pass after final GTR assignment

### Joint Reconstruction

- [ ] Removed in v1 ([intentional change](../port-intentional-changes/ancestral-joint-reconstruction-removed.md))
- [ ] Parsed and exposed in CLI but runtime path is `unimplemented!`

### GTR Model Inference

- [x] GTR inference from data (sparse and dense paths)
- [x] Uninformative root state filtering in dense path ([intentional change](../port-intentional-changes/gtr-uninformative-root-state-filtering.md))
- [x] GTR bootstrapped through temporary JC69 model before replacement
- [ ] Iterative GTR inference (`infer_gtr_iterative()` in v0)

### CLI Options

- [x] `--reconstruct-tip-states` (overwrite ambiguous tips, controls leaf emission)
- [x] `--model` relevant only on marginal paths (parsimony bypasses GTR)
- [ ] `--keep-overhangs` (parsed but not wired)
- [ ] `--zero-based` indexing (parsed but not wired)
- [ ] `--report-ambiguous` (parsed but not wired)
- [ ] `--seed` for reproducibility (parsed but not wired)
- [ ] `--gtr-params` custom GTR parameters (parsed but not wired)
- [ ] `--aa` (parsed but not wired)
- [ ] `--vcf-reference` (parsed, VCF reader not implemented)
- [ ] `--aln` legacy option (parsed but not wired)
- [ ] Tree inference from alignment (help text mentions it, not implemented)

## 2. Clock Inference

### Clock Regression

- [x] Backward regression pass
  - [x] Leaf contribution (date + divergence)
  - [x] Outlier leaf contribution path
  - [x] Internal node aggregation
  - [x] Root assignment
- [x] Forward regression pass for reroot search
- [x] Variance propagation (branch length proportional + offset)
- [x] Covariation mode (`--covariation` for regression)
  - [x] Uses `--sequence-length`
  - [x] Uses `--tip-slack` or default `3.0`
  - [x] Builds custom variance parameters
  - [/] Missing `--sequence-length` does not raise explicit error

### Rerooting

- [x] Node scan for best positive-rate root candidate
- [x] Parent-branch split optimization around best node
- [x] Child-branch split optimization around best node
- [x] Reroot rejects candidates with non-positive inferred clock rate
- [x] `--keep-root` (disable rerooting)

### Reroot Modes (5 modes)

- [x] Least-squares (minimize RTT regression residuals)
- [x] Min-dev (minimize RTT variance)
- [x] Oldest (root at oldest node)
- [x] Clock-filter mode
- [x] MRCA mode

### Edge Split Optimization (3 methods)

- [x] Grid search (configurable grid points)
- [x] Brent's method (configurable tolerance/max-iter)
- [x] Golden section search (configurable tolerance/max-iter)

### Edge Operations for Rerooting

- [x] Snap to source endpoint
- [x] Snap to target endpoint
- [x] Edge split (create new node at split point)
- [x] Invert every edge on path to new root
- [x] Edge merge (remove trivial nodes, remove old root when it becomes trivial)

### Clock Filtering

- [x] Outlier prefilter when `--clock-filter > 0`
- [x] Divergence assignment (root-to-tip accumulation)
- [x] IQD calculation (Q3-Q1)
- [x] Outlier marking (threshold \* IQD)
- [ ] Local clock filter method (`--clock-filter-method=local` in v0)
- [ ] `--prune-outliers` (v0 clock command removes outliers from output tree)

### Clock Model Output

- [x] Rate, intercept, chisq, R-value, hessian, covariance
- [x] Rerooted Newick tree
- [x] Graph JSON and Graphviz DOT
- [x] Clock model JSON
- [x] Regression CSV (populated for all nodes, includes predicted date, clock deviation, outlier state)
- [x] SVG and PNG RTT charts
- [x] Terminal chart printing on TTY

### CLI Options

- [x] `--clock-rate` (fixed rate bypasses estimation)
- [x] `--allow-negative-rate` (for midpoint rooting)
- [x] `--tip-slack` (terminal node excess variance)
- [ ] `--aln` (parsed but not wired)
- [ ] `--vcf-reference` (parsed but not wired)
- [ ] `--gtr` (parsed but not wired)
- [ ] `--gtr-params` (parsed but not wired)
- [ ] `--branch-length-mode` (parsed but not wired)
- [ ] `--method-anc` (parsed but not wired)
- [ ] `--reroot` (parsed but not wired)
- [ ] `--prune-short` (parsed but not wired)
- [ ] `--plot-rtt` (parsed but not wired)
- [ ] `--seed` (parsed but not wired)
- [ ] Tree inference from alignment (help text mentions it, not implemented)

## 3. Timetree Inference

### Initialization

- [x] Load tree from `--tree`
- [x] Optional date table from `--dates`
- [x] Optional alignment when `--branch-length-mode=marginal`
- [x] Alphabet always maps gap to unknown profile (matching v0 `nuc_nogap` behavior)
- [x] Initial node divergences always initialized
- [x] Date assignment fails when no valid dates or fewer than three leaves have valid dates
- [ ] Tree inference from alignment path is `todo!`

### Time Distribution Propagation

- [x] Backward pass (leaf dates to root via convolution)
- [x] Forward pass (root to leaves, refine distributions)
- [x] Bad branch exclusion (outliers, dateless leaves)
- [x] Build branch distributions from partitions when present
- [x] Build point branch distributions from input lengths when partitions absent
- [x] Optional coalescent contribution calculation

### Branch-Length Modes

- [x] `input` (use tree branch lengths, no sequence partitions, point distributions)
- [x] `marginal` (optimize via marginal reconstruction)
  - [x] Dense marginal partition branch
  - [x] Sparse marginal partition branch (compresses sequences first)
  - [x] GTR model selection from `--gtr` flag (named models and inference)
- [ ] `joint` (v0 supports joint optimization)

### Time Marginal Modes

- [x] `never` (joint most-likely times)
- [ ] `always` (marginal every round - parsed but not wired)
- [x] `only-final` (marginal last round for confidence)
  - [x] Final timetree pass after loop
  - [x] Final marginal update when partitions exist
  - [x] Confidence interval extraction from node time distributions
  - [x] Confidence TSV output
- [x] `confidence-only` (v0 variant, implemented as `--confidence` promoting `Never` to `OnlyFinal`)

### Coalescent Models

- [x] Constant Tc (fixed from CLI `--coalescent`)
- [x] Optimized Tc (`--coalescent-opt`, Brent's method, log-space bracket [-20, 2])
  - [x] Re-optimizes constant Tc inside loop for iterations i >= 2
  - [/] No special pre-loop setup unless skyline also active
- [x] Skyline (`--coalescent-skyline`, Nelder-Mead, piecewise linear Tc(t))
  - [x] Pre-loop constant Tc optimization
  - [x] Warning fallback to Tc = 1.0 on failure or non-convergence
  - [x] Final skyline re-optimization after refinement loop
  - [x] Extra final timetree pass unless `--time-marginal=only-final`
- [x] Coalescent contribution per node (survival + merger probability)
- [x] Different multiplication ordering ([intentional change](../port-intentional-changes/coalescent-multiplication-ordering.md))
- [x] Merger rate lambda(t) = k(k-1)/(2\*Tc)
- [x] Branch counting k(t) from node times
- [ ] `--n-branches-posterior` (parsed, returns explicit error - [known issue](../port-known-issues/N-timetree-n-branches-posterior-unimplemented.md))
- [ ] Empirical skyline (v0 `skyline_empirical()` - sliding window without optimization)
- [ ] Skyline confidence intervals (v0 computes via second derivatives)
- [ ] Skyline plot output ([known issue](../port-known-issues/N-timetree-missing-skyline-output.md))

### Polytomy Resolution

- [x] Find polytomy nodes with more than two children
- [x] Greedy pairwise merging (likelihood gain threshold 0.05)
- [x] Per-pair Brent optimization for merge time
- [x] Remove obsolete single-child nodes after resolution
- [x] Reconcile partition topology after tree change
- [/] `--keep-polytomies` (parsed but never read - [known issue](../port-known-issues/N-timetree-dead-cli-flags.md))
- [ ] Stochastic resolution (`--stochastic-resolve` in v0 - [known issue](../port-known-issues/N-timetree-stochastic-polytomy-unimplemented.md))

### Relaxed Clock

- [x] Postorder penalty coefficients (k1, k2)
- [x] Preorder optimal gamma (rate multiplier per branch)
- [x] Coupling parameter (parent-offspring rate correlation)
- [x] Slack parameter (rate deviation penalty)
- [x] Skip relaxed clock when total partition length is zero
- [x] Use first two values as slack and coupling, with defaults when omitted
- [x] `--relax` argument parsing (`num_args = 2` for slack and coupling)
- [ ] Substitution rates output (`substitution_rates.tsv` in v0)

### Refinement Iteration Loop

- [x] Loop stops on convergence or `--max-iter`
- [x] Convergence check uses `n_diff == 0 && n_resolved == 0`
- [x] Skyline mode suppresses early convergence exit
- [x] Relaxed clock application
- [x] Polytomy resolution
- [x] Dirty-tree-aware reconstruction ordering
  - [x] If topology changed, rerun timetree first and marginal second
  - [x] If topology did not change, rerun marginal first and timetree second
- [x] Re-estimate clock model after every iteration without rerooting
- [x] Optional tracelog write on each recorded iteration

### Convergence Tracking

- [x] Sequence change count (n_diff)
- [x] Polytomies resolved count (n_resolved)
- [x] Sequence likelihood (lh_seq)
- [x] Positional likelihood (lh_pos)
- [x] Tracelog CSV output (`--tracelog`)
- [x] Coalescent likelihood in metrics (lh_coal)

### Confidence Intervals

- [x] 90% HPD from marginal posteriors
- [x] `combine_confidence()` quadrature sum
- [x] `Distribution::quantile()` inverse CDF
- [x] Rate susceptibility analysis (`compute_rate_susceptibility()`)
- [x] `date_uncertainty_due_to_rate()` via probit function (erfinv)
- [x] Rate susceptibility activated via `--confidence` with `--covariation` or `--clock-std-dev`
- [x] `--confidence` promotion of `time_marginal` from `never` to `only-final`
- [x] `--covariation` wired into timetree clock regression

### Output

- [x] Timetree Newick
- [x] Timetree Nexus
- [x] Clock model JSON with `timetree.*` basename
- [x] Confidence TSV
- [ ] Node dates TSV (`write_node_dates()` is `todo!()` - [known issue](../port-known-issues/N-timetree-node-dates-output-unimplemented.md))
- [ ] Substitution rates TSV (v0 writes `substitution_rates.tsv` when `--relax` is used)
- [ ] Auspice JSON (v0 writes `auspice_tree.json`)
- [ ] Outliers TSV (v0 writes `outliers.tsv`)
- [ ] Tracelog run (v0 `tracelog_run()` with detailed per-iteration state)
- [ ] Plotting (`--plot-tree`, `--plot-rtt` - parsed, return explicit error - [known issue](../port-known-issues/N-timetree-plot-unimplemented.md))

### CLI Options (Parsed but Not Wired)

- [x] `--clock-std-dev` (rate susceptibility)
- [x] `--confidence` (time_marginal promotion, rate susceptibility)
- [x] `--covariation` (GLS clock regression, rate susceptibility)
- [x] `--tip-slack` (covariation variance computation)
- [ ] `--n-iqd`
- [ ] `--reroot`
- [ ] `--tip-labels` / `--no-tip-labels`
- [x] `--gtr` (model selection: named models and inference)
- [ ] `--gtr-params` (parsed but not wired)
- [ ] `--method-anc`
- [ ] `--aa`
- [ ] `--keep-overhangs`
- [ ] `--zero-based`
- [ ] `--reconstruct-tip-states`
- [ ] `--report-ambiguous`
- [ ] `--seed`
- [ ] `--gen-per-year` (generations per year for N_e estimation)
- [ ] `--aln` legacy option

## 4. Homoplasy Analysis

- [ ] **Entire command unimplemented in v1** (runtime path is `unimplemented!`)
- [x] CLI shape exists, inherits all `ancestral` arguments through flattening

### v0 Features (All Missing)

- [ ] Mutation multiplicity distribution
- [ ] Position hit count distribution
- [ ] Poisson comparison (expected vs observed)
- [ ] Top-N homoplasic mutations display
- [ ] `--detailed` (terminal branch mutations, strains with homoplasies)
- [ ] `--drms` (drug resistance mutation annotation)
- [ ] `--const` (constant sites correction)
- [ ] `--rescale` (branch length rescaling)
- [ ] `-n` (number of mutations/nodes to display)

## 5. Mugration (Discrete Trait Reconstruction)

- [x] **Implemented** (marginal reconstruction complete, ~80%)

### Input Processing

- [x] Tree loading from `--tree`
- [x] Discrete state table from `--states`
- [x] Attribute column selection with `--attribute`
- [x] Custom taxon-name column with `--name-column`
- [x] Auto-detect taxon-name column from `name`, `strain`, or `accession`
- [x] Weights file parsing and validation
- [x] Merge categories observed in states and weights files
- [x] Warn on categories missing from weights file
- [x] Enforce `--missing-weights-threshold`
- [x] Fill missing weights with mean observed weight
- [x] Normalize weights to sum to one
- [x] Build discrete attribute alphabet excluding missing-data token
- [x] Add synthetic missing-data symbol after alphabet build
- [x] Reject datasets with fewer than two non-missing states

### Reconstruction

- [x] GTR model construction (uniform rates, optional weight-based equilibrium frequencies)
- [x] Discrete marginal reconstruction (backward pass, forward pass)
- [x] Missing data handling (uniform prior for tips with `"?"` traits)
- [x] Trait assignment from argmax of posterior profiles
- [x] Confidence profile extraction (`get_confidence()`)

### Output

- [x] `traits.csv` (per-node trait assignments, all nodes)
- [x] `annotated_tree.nexus` (Newick with trait annotations)
- [x] `annotated_tree.nwk` (Newick with NHX-style annotations)
- [x] `gtr.json` (GTR model parameters)

### Additional Features

- [x] Iterative GTR inference ([golden master parity tracked](../port-known-issues/M-mugration-iterative-gtr.md))
- [x] Sampling bias correction (`--sampling-bias-correction`)
- [x] Confidence CSV output (`--confidence`)
- [x] `--pc` pseudo-counts

## 6. ARG Inference

- [ ] **Not present in v1**

### v0 Features (All Missing)

- [ ] Two-tree input (`--trees`, `--alignments`)
- [ ] MCC (most compatible clades) file input
- [ ] Per-tree timetree inference
- [ ] Combined ARG output

## 7. Branch Length Optimization

v0 has no standalone optimize command. `TreeAnc.optimize_tree()` and
`optimize_tree_marginal()` perform branch length optimization inline during
timetree and ancestral workflows. v1 extracts this into a standalone `optimize`
command and a reusable `PartitionOptimizeOps` trait system shared with timetree.

### Standalone Command (v1-only)

- [x] `optimize` command with `--tree`, `--aln`, `--outdir`
- [x] Mixed dense and sparse partition setup (one of each, from same alignment)
- [x] Initial branch-length guess from observed mutation counts (excludes deletion and ambiguous positions)
- [x] `--branch-length-initial-guess` flag: `auto` (selective fill, treats zero BL as invalid when indels present), `always` (overwrite all), `never` (error on missing)
- [x] Output annotated Newick and Nexus trees
- [x] GTR parameters written to JSON

### Per-Edge Likelihood (v0 parity via different method)

v0 uses Brent's method (`scipy.optimize.minimize_scalar`) in sqrt(t) space with
Hamming distance bracket. v1 uses Newton's method with analytical derivatives
in t-space, falling back to grid search. Both use eigendecomposition-based
likelihood (`expQt = V diag(exp(lambda*t)) V_inv`).

- [x] Eigenvalue-space coefficient caching (dense: `msg.dot(V) * msg.dot(V_inv.T)`, sparse: per-site with multiplicity)
- [x] Analytical first and second derivatives (v1-only, v0 uses derivative-free Brent)
- [x] Newton's method with clamped step (max 10 inner iterations, step in `[-1.0, bl]`, absolute tolerance floor 1e-8 subs/site)
- [x] Grid search fallback when second derivative >= 0 (100 points, log-spaced grid with 0.5 subs/site minimum upper bound)
- [x] Zero branch length short-circuit (combined likelihood > 0.01 and derivative < 0 at zero)
- [x] `compute_derivatives` flag to skip derivative computation for log-likelihood-only evaluation
- [x] Collect dense contribution (`PartitionMarginalDense`)
- [x] Collect sparse contribution (`PartitionMarginalSparse`, multiplicity-weighted)
- [x] Unified mixed-partition evaluation (`evaluate_mixed()` sums metrics across partition types)
- [x] sqrt(t) reparameterization (Brent in s=sqrt(t) space, matching v0; also ln(t) space variant)
- [ ] Regularization penalty for profile-based optimization (v0: `exp(t^4/10000)` prevents unbounded growth)
- [ ] Hamming distance fallback when optimization fails (v0 falls back to observed distance)

### Convergence Loop

- [x] Iterative marginal reconstruction + optimization loop bounded by `--max-iter`
- [x] Early stop when absolute likelihood change is below `--dp`
- [ ] Damping in marginal loop (v0: `new*(1-d^i) + old*d^i` with d=0.75)
- [ ] Progressive per-iteration tolerance tightening (v0: `tol = 1e-8 + 0.01^(i+1)`, coarse early, tight late)
- [ ] Bifurcating root special handling (v0 optimizes combined root-children length, preserves ratio)
- [ ] Convergence by sequence change count (v0 joint mode: stops when zero nucleotides change)

### GTR Integration

- [x] `--model` flag wired through `get_gtr_sparse()`/`get_gtr_dense()` for all named models and inference
- [ ] GTR inference integrated into optimization loop (v0: `infer_gtr` parameter re-estimates model per iteration)
- [ ] GTR rate optimization (v0: `optimize_gtr_rate()` optimizes overall substitution rate mu)

### Reuse by Timetree

- [x] `PartitionOptimizeOps` trait implemented by both dense and sparse partitions
- [x] `collect_edge_contributions()` gathers contributions for one edge across partition types
- [x] `compute_branch_length_distribution()` evaluates log-likelihood on grid, converts to time-domain distribution
- [x] `evaluate_mixed_log_lh_only()` for grid evaluation without derivatives

### Branch Length Modes (v0 has 3 modes, v1 implements marginal only)

- [x] Marginal mode (profile-based optimization with full probability vectors)
- [ ] Joint mode (removed in v1, see [intentional change](../port-intentional-changes/ancestral-joint-reconstruction-removed.md))
- [ ] Input mode (no optimization, Poisson/Gaussian approximation for short/long branches)
- [ ] Auto-detection of branch length mode (v0: defaults to `input` if max branch > 0.1, else `joint`)

### Short Branch Handling

- [x] Zero-branch-length derivative check before optimization
- [ ] Short branch pruning after optimization (v0 joint mode: inside loop; v0 marginal mode: after loop; prunes when `bl < 0.1 * one_mutation` and P(zero) > 0.1)
- [ ] MIN_BRANCH_LENGTH floor for GTR calculations (v0: `1e-3 * one_mutation`)

### Current Limitations

- [/] Command always builds one sparse and one dense partition from the same full alignment
- [ ] `--dense` (parsed but not wired, `infer_dense()` is a stub returning false)
- [ ] Separate dense-only and sparse-only command modes not exposed
- [ ] Standalone `run_optimize_sparse()` zero-branch inconsistency

## 8. Pruning (v1-Only Command)

- [x] **Prune command** (no v0 standalone equivalent)
- [x] Tree input from `--tree`
- [x] Optional alignment input from `--aln`
- [x] Output pruned Newick and Nexus trees

### Pruning Criteria

- [x] Short branch pruning (`--prune-short` threshold)
- [x] Empty branch pruning (`--prune-empty`, requires alignment)
- [x] Name-based pruning (`--prune-nodes-list`, `--prune-nodes-list-file`)
- [x] Custom delimiters for node lists

### Validation

- [x] `--prune-empty` without `--aln` returns explicit error
- [x] `--prune-empty` loads sparse partitions and compresses sequences
- [x] Mutation counts come from sparse edge substitutions

### Pruning Workflow

- [x] Two-pass strategy
  - [x] Internal-edge collapse pass
  - [x] Graph rebuild
  - [x] Leaf-removal pass
  - [x] Graph rebuild
- [x] Internal edge collapsed when: shorter than threshold, empty of substitutions, or target node name selected
- [x] Leaf pass removes only selected leaf names
- [x] Recursive parent collapse (childless parent removal)
- [x] Edge collapse with data merging (branch length sum, mutation union)

## 9. GTR Substitution Models

### Implemented Models

- [x] Named nucleotide models: JC69, K80, F81, HKY85, T92, TN93
- [x] Amino acid model: JTT92

### Model Inference

- [x] Model inference from data (sparse and dense paths)
- [x] Eigendecomposition (symmetric, caching)
- [x] expQt (matrix exponential for transition probabilities)
- [x] Equilibrium frequencies (pi vector, normalized)
- [x] Exchangeability matrix (W, symmetric, normalized by avg_transition)
- [x] Rate scaling (mu parameter)
- [x] JSON output (GtrOutput struct with model type, parameters)

### Not Implemented

- [ ] Custom GTR from file (`--custom-gtr` / `GTR.from_file()` in v0)
- [ ] Random GTR generation (`GTR.random()` in v0, used for testing)
- [ ] GTR save to file (`save_to_npz()` in v0)
- [ ] Site-specific models (multi-dimensional pi/W in v0)
- [ ] Branch length optimization at GTR level (v0 `optimal_t()` - moved to partition layer in v1)
- [ ] Sequence probability at GTR level (v0 `prob_t()` - moved to partition layer in v1)
- [ ] State pair compression (v0 `state_pair()` - moved to partition layer in v1)

## 10. Alphabet System

- [x] Nucleotide alphabet (A, C, G, T with gap handling)
- [x] Amino acid alphabet (20 AAs, with/without stop codon)
- [x] Gap always mapped to unknown profile (no toggle needed for nogap alphabets)
- [x] IUPAC ambiguity codes (nucleotide: R, Y, S, W, K, M, D, H, B, V, N, X)
- [x] Amino acid ambiguity (X, B, Z)
- [x] Profile maps (ambiguity code to probability vector, in primitives layer)
- [ ] Alphabet guessing from data (v0 `guess_alphabet()`: >90% nuc threshold)

## 11. Probability Distributions

### Implemented

- [x] Analytical Gaussian (PDF, product, convolution)
- [x] Analytical Exponential (PDF, convolution, special case a=b)
- [x] Gaussian-Exponential convolution
- [x] Grid distributions (uniform grid, linear interpolation)
- [x] ScaledArray pattern (normalized values + log-scale factor)

### Not Implemented

- [ ] Unified Distribution class (v0: wraps scipy.interpolate.interp1d)
- [ ] Delta functions (point masses, v0 `Distribution.delta_function()`)
- [ ] Distribution multiplication (v0 `Distribution.multiply()`)
- [ ] Distribution division (v0 `Distribution.divide()`)
- [ ] Numerical integration (v0 Simpson's rule, trapezoidal)
- [ ] FFT transform (v0 `Distribution.fft()`)
- [ ] FWHM calculation (full width half maximum)
- [ ] Effective support (threshold-based range)
- [ ] Grid refinement (v0 `_adjust_grid()` adaptive)
- [ ] BranchLenInterpolator (v0: node-specific branch length probability)
  - [ ] Input mode (Gaussian/Poisson approximation)
  - [ ] Marginal mode (from profile pairs)
  - [ ] Joint mode (from compressed state pairs)
  - [ ] Gamma rescaling (for relaxed clock)
  - [ ] Adaptive grid construction (log near zero, quadratic tails)

## 12. Representation and Partition System

### Partitions

- [x] PartitionFitch (Fitch parsimony, sparse storage)
- [x] PartitionMarginalDense (full probability vectors at all positions)
- [x] PartitionMarginalSparse (variable positions only)
- [x] Partition traits (PartitionMarginal, PartitionMarginalOps, PartitionCompressed, HasLogLh)

### Payloads

- [x] Dense payloads (DenseNodePartition, DenseEdgePartition, DenseSeqDis)
- [x] Sparse payloads (SparseNodePartition, SparseEdgePartition, MarginalSparseSeqDistribution)
- [x] Node payloads (NodeAncestral, NodeTimetree with time distributions)
- [x] Edge payloads (EdgeAncestral, EdgeTimetree with branch distributions and messages)

### Operations

- [x] Mutation composition (edge merge: A5G + G5T = A5T with cancellation)
- [x] Edge inversion (swap ref/qry, reverse messages for rerooting)
- [x] Reroot operations on partitions (split/merge edges, update messages)

### Architecture

- [x] Dense/sparse architecture ([intentional change](../port-intentional-changes/sequence-representation-dense-sparse.md))
- [x] Partition system ([intentional change](../port-intentional-changes/partition-system-architecture.md))
- [x] Graph-based phylogenetic representation ([intentional change](../port-intentional-changes/graph-based-phylogenetic-representation.md))

## 13. Sequence Primitives

- [x] AsciiChar (validated u8, type-safe character)
- [x] Seq (Vec\<AsciiChar\> with string-like operations, concatenation, repetition)
- [x] BitSet128/StateSet (u128 bitfield for character sets, set operations)
- [x] Sub (substitution mutation: pos + ref + qry, invertible, parseable)
- [x] InDel (insertion/deletion: range + seq + direction, invertible)
- [x] Composition (character frequency tracking with mutation updates)
- [x] find_char_ranges (contiguous range detection for gaps, unknowns)
- [x] compute_divs (root-to-tip divergence calculation)

## 14. I/O Formats

### Input Formats

- [x] FASTA (plain text)
- [x] Compressed FASTA (gzip, xz, bzip2, zstd)
- [x] Multiple FASTA files (concatenated)
- [x] Stdin input
- [x] Newick tree
- [x] Nexus tree
- [x] Phylip tree
- [x] CSV/TSV dates
- [x] CSV/TSV discrete states (mugration)
- [ ] VCF input (variant call format)
- [ ] Compressed VCF (.vcf.gz)
- [ ] Custom GTR model file

### Output Formats

- [x] FASTA (ancestral sequences)
- [x] Newick (annotated tree)
- [x] Nexus (annotated tree)
- [x] JSON (GTR model, clock model, graph)
- [x] CSV (clock regression, confidence intervals)
- [x] SVG/PNG charts (clock regression)
- [x] Graphviz DOT
- [ ] VCF output (v0 writes .vcf for VCF inputs)
- [ ] Auspice JSON (Nextstrain visualization)
- [ ] Skyline TSV/plot
- [ ] Substitution rates TSV
- [ ] Outliers TSV
- [ ] Branch mutations TXT
- [ ] Sequence evolution model TXT

### v1-Only Formats

- [x] PhyloXML
- [x] Usher MAT (mutation-annotated tree)
- [x] YAML serialization
- [x] Compressed FASTA output
- [x] Streaming readers/writers with automatic decompression

## 15. Date Parsing

- [x] Float dates (e.g. 2012.15)
- [x] ISO format (YYYY-MM-DD)
- [x] Date ranges ([2013.2:2013.7])
- [x] Custom name/date column selection
- [x] Imprecise dates (YYYY-XX-XX) ([known issue: upper bound not capped](../port-known-issues/N-dates-imprecise-upper-bound-not-capped.md))
- [/] Auto-detect column names ([known issue: detection gaps](../port-known-issues/M-dates-column-auto-detection-gaps.md))

## Statistics

Counts represent top-level feature groups (section headings and their direct children), not individual leaf checkboxes. Nested sub-items are rolled up into their parent feature.

| Domain         | v0 Features | v1 Implemented | v1 Missing | v1 Only |
| -------------- | ----------- | -------------- | ---------- | ------- |
| ancestral      | 16          | 10             | 6          | 0       |
| clock          | 17          | 14             | 3          | 0       |
| timetree       | 30          | 21             | 9          | 0       |
| homoplasy      | 9           | 0              | 9          | 0       |
| mugration      | 10          | 10             | 0          | 0       |
| arg            | 4           | 0              | 4          | 0       |
| optimize       | 0           | 0              | 0          | 8       |
| prune          | 0           | 0              | 0          | 7       |
| gtr            | 12          | 7              | 5          | 0       |
| alphabet       | 6           | 5              | 1          | 0       |
| distribution   | 15          | 5              | 10         | 0       |
| representation | 12          | 12             | 0          | 0       |
| primitives     | 7           | 7              | 0          | 0       |
| io-input       | 11          | 9              | 2          | 0       |
| io-output      | 13          | 6              | 7          | 4       |
| dates          | 5           | 4              | 1          | 0       |

## Cross-Command Notes

- `ancestral`, `timetree`, `homoplasy` expose `MethodAncestral`
- `clock` and `timetree` share the same clock-regression and reroot machinery
- `timetree` reuses `ancestral::marginal` and `clock::*` internals rather than maintaining a separate reconstruction stack
- `optimize` and `prune` are v1-only commands
- Several help strings describe planned behavior not wired in current command code

## Cross-References

- [Algorithm Inventory](../port-algo-inventory/_index.md) - implementation details
- [Intentional Changes](../port-intentional-changes/_index.md) - deliberate deviations
- [Known Issues](../port-known-issues/_index.md) - bugs and missing features
- [Test Inventory](../port-test-inventory/_index.md) - test coverage
