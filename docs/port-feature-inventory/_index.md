# Feature Inventory - TreeTime v0/v1

> **Status**: reviewed against source code
>
> **Valid for commit**: `404d8b77` on 2026-03-10
>
> **Scope**: v0/v1 feature parity tracking. Combines domain taxonomy with CLI wiring status.

## Legend

- `[x]` implemented in v1 with v0 parity (or v1-only feature)
- `[~]` partial implementation, stubbed, or only partly wired
- `[ ]` missing, parsed but not wired, or documented but unimplemented

## Command Map

| Command     | Status | Notes                                 |
| ----------- | ------ | ------------------------------------- |
| `ancestral` | [x]    | Parsimony and marginal reconstruction |
| `clock`     | [x]    | Regression and rerooting              |
| `timetree`  | [x]    | Full inference pipeline               |
| `optimize`  | [x]    | v1-only branch-length optimization    |
| `prune`     | [x]    | v1-only tree pruning                  |
| `mugration` | [~]    | Preprocessing only, no reconstruction |
| `homoplasy` | [ ]    | Unimplemented                         |

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
- [~] GTR bootstrapped through temporary JC69 model before replacement
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
  - [~] Missing `--sequence-length` does not raise explicit error

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
- [x] Alphabet setup changes `treat_gap_as_unknown` based on alignment presence and dense mode
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
  - [~] GTR hardcoded to JC69 for timetree command
- [ ] `joint` (v0 supports joint optimization)

### Time Marginal Modes

- [x] `never` (joint most-likely times)
- [x] `always` (marginal every round)
- [x] `only-final` (marginal last round for confidence)
  - [x] Final timetree pass after loop
  - [x] Final marginal update when partitions exist
  - [x] Confidence interval extraction from node time distributions
  - [x] Confidence TSV output
- [~] `never` and `always` parsed as separate enum values but only `only-final` has distinct behavior
- [ ] `confidence-only` (v0 variant)

### Coalescent Models

- [x] Constant Tc (fixed from CLI `--coalescent`)
- [x] Optimized Tc (`--coalescent-opt`, Brent's method, log-space bracket [-20, 2])
  - [x] Re-optimizes constant Tc inside loop for iterations i >= 2
  - [~] No special pre-loop setup unless skyline also active
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
- [x] `--keep-polytomies` (disable resolution)
- [ ] Stochastic resolution (`--stochastic-resolve` in v0 - [known issue](../port-known-issues/N-timetree-stochastic-polytomy-unimplemented.md))

### Relaxed Clock

- [x] Postorder penalty coefficients (k1, k2)
- [x] Preorder optimal gamma (rate multiplier per branch)
- [x] Coupling parameter (parent-offspring rate correlation)
- [x] Slack parameter (rate deviation penalty)
- [x] Skip relaxed clock when total partition length is zero
- [x] Use first two values as slack and coupling, with defaults when omitted
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
- [ ] Coalescent likelihood in metrics ([known issue](../port-known-issues/M-timetree-coalescent-likelihood-stub.md))

### Confidence Intervals

- [x] 95% from marginal posteriors
- [ ] Rate susceptibility analysis (TODO in code)
- [ ] `--vary-rate` sensitivity (parsed, returns explicit error - [known issue](../port-known-issues/H-timetree-vary-rate-unimplemented.md))

### Output

- [x] Timetree Newick
- [x] Timetree Nexus
- [x] Clock model JSON with `timetree.*` basename
- [x] Confidence TSV
- [ ] Auspice JSON (v0 writes `auspice_tree.json`)
- [ ] Outliers TSV (v0 writes `outliers.tsv`)
- [ ] Tracelog run (v0 `tracelog_run()` with detailed per-iteration state)
- [ ] Plotting (`--plot-tree`, `--plot-rtt` - parsed, return explicit error - [known issue](../port-known-issues/N-timetree-plot-unimplemented.md))

### CLI Options (Parsed but Not Wired)

- [ ] `--clock-std-dev`
- [ ] `--confidence`
- [ ] `--n-iqd`
- [ ] `--reroot`
- [ ] `--tip-labels` / `--no-tip-labels`
- [ ] `--covariation`
- [ ] `--gtr` / `--gtr-params`
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

- [~] **Partially implemented** (preprocessing only, ~30%)

### Implemented

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

### Not Implemented

- [ ] Tree loading or reconstruction work
- [ ] GTR model construction from discrete traits (commented out)
- [ ] Ancestral state reconstruction on discrete traits
- [ ] Sampling bias correction (`--sampling-bias-correction`)
- [ ] Confidence output (`--confidence`)
- [ ] Output files (GTR.txt, annotated_tree.nexus, confidence.csv)

### CLI Options (Parsed but Not Wired)

- [ ] `--tree`
- [ ] `--confidence`
- [ ] `--pc`
- [ ] `--sampling-bias-correction`
- [ ] `--outdir`
- [ ] `--seed`

## 6. ARG Inference

- [ ] **Not present in v1**

### v0 Features (All Missing)

- [ ] Two-tree input (`--trees`, `--alignments`)
- [ ] MCC (most compatible clades) file input
- [ ] Per-tree timetree inference
- [ ] Combined ARG output

## 7. Branch Length Optimization (v1-Only Command)

- [x] **Optimize command** (no v0 standalone equivalent)
- [x] Tree input from `--tree`
- [x] Alignment input from `--aln`
- [x] Mixed dense and sparse partition setup
- [x] Initial branch-length guess from observed differences
- [x] Output annotated Newick and Nexus trees

### Per-Edge Optimization

- [x] Collect dense contribution
- [x] Collect sparse contribution
- [x] Newton's method (when second derivative < 0)
- [x] Grid search fallback (when Newton unavailable)
- [x] Zero branch length optimum short-circuit
- [x] Clamp Newton step to avoid overshoot
- [x] Coefficient extraction (dense: eigenvector products, sparse: variable position collection)

### Convergence

- [x] Iterative likelihood loop bounded by `--max-iter`
- [x] Early stop when absolute likelihood change is below `--dp`

### Current Limitations

- [~] Command always builds one sparse and one dense partition from the same full alignment
- [~] GTR hardcoded to JC69
- [ ] Separate dense-only and sparse-only command modes not exposed
- [ ] `--model` (parsed but not wired)
- [ ] `--dense` (parsed but not wired)

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
- [x] Gap-as-unknown toggle (`treat_gap_as_unknown` per partition)
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
- [x] Graph-based tree structure ([intentional change](../port-intentional-changes/tree-structure-graph-based.md))

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
- [ ] Imprecise dates (YYYY-XX-XX)
- [ ] Auto-detect column names (v0 tries `name`, `strain`, `accession`)

## Statistics

| Domain         | v0 Features | v1 Implemented | v1 Missing | v1 Only |
| -------------- | ----------- | -------------- | ---------- | ------- |
| ancestral      | 16          | 10             | 6          | 0       |
| clock          | 17          | 14             | 3          | 0       |
| timetree       | 30          | 20             | 10         | 0       |
| homoplasy      | 9           | 0              | 9          | 0       |
| mugration      | 10          | 4              | 6          | 0       |
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
| dates          | 5           | 3              | 2          | 0       |

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
