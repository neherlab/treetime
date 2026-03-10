# Feature Inventory - TreeTime v0/v1

> **Status**: In progress
>
> AI-generated via `feature-inventory` skill at commit `a763ad09` on 2026-03-10. Verify against source code.
>
> `[x]` = implemented in v1 with v0 parity. `[ ]` = missing or partial in v1. `[v1]` = v1-only feature (no v0 equivalent).

## Feature Taxonomy

### 1. Ancestral Reconstruction

- [x] **Fitch parsimony** (sparse)
  - [x] Backward pass (intersection/union of child state sets)
  - [x] Forward pass (resolve ambiguities top-down)
  - [x] Indel tracking (insertions, deletions)
  - [x] Composition tracking (character counts)
  - [x] Deterministic root state resolution ([intentional change](../port-intentional-changes/ancestral-fitch-deterministic-root-state.md))
- [x] **Marginal reconstruction** (dense)
  - [x] Backward pass (log-space message multiplication)
  - [x] Forward pass (parent-child message propagation)
  - [x] Root equilibrium frequency weighting
  - [x] Log-space normalization (logsumexp, matches v0)
  - [x] expQt matrix propagation
  - [ ] Profile sampling (`sample_from_prof` - v0 supports stochastic sampling, v1 argmax only)
- [x] **Marginal reconstruction** (sparse)
  - [x] Backward pass (variable positions only)
  - [x] Forward pass (variable positions only)
  - [x] Fixed-position per-character profiles
  - [x] Mutation composition from edge substitutions
  - [x] Message combine with variability threshold (EPS=1e-4)
- [ ] **Joint reconstruction** - removed in v1 ([intentional change](../port-intentional-changes/ancestral-joint-reconstruction-removed.md))
- [x] **GTR model inference** from data (sparse and dense paths, uninformative root state filtering in dense path - [intentional change](../port-intentional-changes/gtr-uninformative-root-state-filtering.md))
- [ ] **Iterative GTR inference** (`infer_gtr_iterative()` in v0)
- [x] **Sequence reconstruction** (depth-first traversal, parent + mutations)
- [x] `--reconstruct-tip-states` (overwrite ambiguous tips)
- [ ] `--keep-overhangs` (parsed but not wired up in v1)
- [ ] `--zero-based` indexing (parsed but not wired up in v1)
- [ ] `--report-ambiguous` (parsed but not wired up in v1)
- [ ] `--seed` for reproducibility (parsed but not wired up in v1)
- [ ] `--gtr-params` custom GTR parameters (parsed but not used in v1)
- [ ] VCF input (`--vcf-reference` parsed but VCF reader not implemented)
- [ ] Tree inference from alignment (fasttree/iqtree/raxml integration)

### 2. Clock Inference

- [x] **Clock regression** (backward/forward message passing with ClockSet)
  - [x] Leaf contribution (date + divergence)
  - [x] Internal node aggregation
  - [x] Root handling
  - [x] Variance propagation (branch length proportional + offset)
- [x] **Rerooting** (5 modes)
  - [x] Least-squares (minimize RTT regression residuals)
  - [x] Min-dev (minimize RTT variance)
  - [x] Oldest (root at oldest node)
  - [x] Clock-filter mode
  - [x] MRCA mode
  - [x] `--keep-root` (disable rerooting)
- [x] **Edge split optimization** (3 methods)
  - [x] Grid search (configurable grid points)
  - [x] Brent's method (configurable tolerance/max-iter)
  - [x] Golden section search (configurable tolerance/max-iter)
- [x] **Edge operations for rerooting**
  - [x] Edge split (create new node at split point)
  - [x] Edge merge (remove trivial nodes)
  - [x] Snap-to-endpoint logic (split near 0.0 or 1.0)
- [x] **Clock filtering** (IQD-based outlier detection)
  - [x] Divergence assignment (root-to-tip accumulation)
  - [x] IQD calculation (Q3-Q1)
  - [x] Outlier marking (threshold \* IQD)
  - [ ] Local clock filter method (`--clock-filter-method=local` in v0)
- [x] **Clock model output** (rate, intercept, chisq, R-value, hessian, covariance)
- [x] **Fixed clock rate** (`--clock-rate` bypasses estimation)
- [x] **Covariation mode** (`--covariation` for regression)
- [x] `--allow-negative-rate` (for midpoint rooting)
- [x] `--tip-slack` (terminal node excess variance)
- [ ] `--prune-outliers` (v0 clock command removes outliers from output tree)
- [ ] Tree inference from alignment
- [x] **Output**: rerooted tree, clock model JSON, regression CSV, SVG/PNG charts

### 3. Timetree Inference

- [x] **Time distribution propagation**
  - [x] Backward pass (leaf dates to root via convolution)
  - [x] Forward pass (root to leaves, refine distributions)
  - [x] Bad branch exclusion (outliers, dateless leaves)
  - [x] Branch distribution creation (marginal mode: grid-based)
  - [x] Branch distribution creation (input mode: point distributions)
- [x] **Time marginal modes**
  - [x] `never` (joint most-likely times)
  - [x] `always` (marginal every round)
  - [x] `only-final` (marginal last round for confidence)
  - [ ] `confidence-only` (v0 variant)
- [x] **Coalescent models**
  - [x] Constant Tc (fixed from CLI)
  - [x] Optimized Tc (Brent's method, log-space bracket [-20, 2])
  - [x] Skyline (Nelder-Mead, piecewise linear Tc(t))
  - [x] Coalescent contribution per node (survival + merger probability, different multiplication ordering - [intentional change](../port-intentional-changes/coalescent-multiplication-ordering.md))
  - [x] Merger rate lambda(t) = k(k-1)/(2\*Tc)
  - [x] Branch counting k(t) from node times
  - [ ] Posterior branch counting (`--n-branches-posterior`, [known issue](../port-known-issues/N-timetree-n-branches-posterior-unimplemented.md))
  - [ ] Empirical skyline (v0 `skyline_empirical()` - sliding window without optimization)
  - [ ] Skyline confidence intervals (v0 computes via second derivatives)
  - [ ] Skyline plot output ([known issue](../port-known-issues/N-timetree-missing-skyline-output.md))
- [x] **Polytomy resolution**
  - [x] Greedy pairwise merging (likelihood gain threshold 0.05)
  - [x] Brent optimization for merge time
  - [x] Post-resolution cleanup (remove single-child nodes)
  - [ ] Stochastic resolution (`--stochastic-resolve` in v0, [known issue](../port-known-issues/N-timetree-stochastic-polytomy-unimplemented.md))
  - [x] `--keep-polytomies` (disable resolution)
- [x] **Relaxed clock**
  - [x] Postorder penalty coefficients (k1, k2)
  - [x] Preorder optimal gamma (rate multiplier per branch)
  - [x] Coupling parameter (parent-offspring rate correlation)
  - [x] Slack parameter (rate deviation penalty)
  - [ ] Substitution rates output (`substitution_rates.tsv` in v0)
- [x] **Convergence tracking**
  - [x] Sequence change count (n_diff)
  - [x] Polytomies resolved count (n_resolved)
  - [x] Sequence likelihood (lh_seq)
  - [x] Positional likelihood (lh_pos)
  - [ ] Coalescent likelihood in metrics ([known issue](../port-known-issues/M-timetree-coalescent-likelihood-stub.md))
  - [x] Tracelog CSV output (`--tracelog`)
- [x] **Confidence intervals** (95% from marginal posteriors)
  - [ ] Rate susceptibility analysis (TODO in code)
  - [ ] `--vary-rate` sensitivity ([known issue](../port-known-issues/H-timetree-vary-rate-unimplemented.md))
- [x] **Refinement iteration loop**
  - [x] Relaxed clock application
  - [x] Polytomy resolution
  - [x] Dirty-tree-aware reconstruction ordering
  - [x] Clock model re-estimation (without rerooting)
- [ ] `--branch-length-mode=joint` (v0 supports joint optimization)
- [x] `--branch-length-mode=input` (use tree branch lengths)
- [x] `--branch-length-mode=marginal` (optimize via marginal reconstruction)
- [x] `--max-iter` (iteration limit)
- [ ] `--gen-per-year` (generations per year for N_e estimation, v0 only)
- [ ] **Plotting** (`--plot-tree`, `--plot-rtt` - [known issue](../port-known-issues/N-timetree-plot-unimplemented.md))
- [ ] **Auspice JSON output** (v0 writes `auspice_tree.json`)
- [ ] **Outliers TSV output** (v0 writes `outliers.tsv`)
- [ ] **Tracelog run** (v0 `tracelog_run()` with detailed per-iteration state)
- [x] **Output**: Newick, Nexus, clock model JSON, confidence TSV

### 4. Homoplasy Analysis

- [ ] **Entire command unimplemented in v1**
- v0 features:
  - [ ] Mutation multiplicity distribution
  - [ ] Position hit count distribution
  - [ ] Poisson comparison (expected vs observed)
  - [ ] Top-N homoplasic mutations display
  - [ ] `--detailed` (terminal branch mutations, strains with homoplasies)
  - [ ] `--drms` (drug resistance mutation annotation)
  - [ ] `--const` (constant sites correction)
  - [ ] `--rescale` (branch length rescaling)
  - [ ] `-n` (number of mutations/nodes to display)

### 5. Mugration (Discrete Trait Reconstruction)

- [ ] **Partially implemented in v1** (~30%)
- [x] Weights file parsing and validation
- [x] Missing value handling (mean weight for absent values)
- [x] Alphabet construction from discrete attributes
- [x] Missing data string filtering
- [ ] GTR model construction from discrete traits (commented out)
- [ ] Ancestral state reconstruction on discrete traits
- [ ] Sampling bias correction (`--sampling-bias-correction`)
- [ ] Confidence output (`--confidence`)
- [ ] Auto-detect column names (v0 tries `name`, `strain`, `accession`)
- [ ] Output: GTR.txt, annotated_tree.nexus, confidence.csv

### 6. ARG Inference

- [ ] **Not present in v1**
- v0 features:
  - [ ] Two-tree input (`--trees`, `--alignments`)
  - [ ] MCC (most compatible clades) file input
  - [ ] Per-tree timetree inference
  - [ ] Combined ARG output

### 7. Branch Length Optimization (v1-only command)

- [v1] **Optimize command** (no v0 standalone equivalent)
  - [v1] Dense partition optimization
  - [v1] Sparse partition optimization
  - [v1] Mixed (dense + sparse) optimization
  - [v1] Newton's method (when second derivative < 0)
  - [v1] Grid search fallback (when Newton unavailable)
  - [v1] Zero branch length check (likelihood + derivative test)
  - [v1] Coefficient extraction (dense: eigenvector products, sparse: variable position collection)
  - [v1] Convergence loop with likelihood threshold (`--dp`)
  - [v1] `--max-iter` iteration limit

### 8. Pruning (v1-only command)

- [v1] **Prune command** (no v0 standalone equivalent)
  - [v1] Short branch pruning (`--prune-short` threshold)
  - [v1] Empty branch pruning (`--prune-empty`, requires alignment)
  - [v1] Name-based pruning (`--prune-nodes-list`, `--prune-nodes-list-file`)
  - [v1] Custom delimiters for node lists
  - [v1] Two-pass strategy (internal nodes first, then leaves)
  - [v1] Recursive parent collapse (childless parent removal)
  - [v1] Edge collapse with data merging (branch length sum, mutation union)

### 9. GTR Substitution Models

- [x] **Named nucleotide models**: JC69, K80, F81, HKY85, T92, TN93
- [x] **Amino acid model**: JTT92
- [x] **Model inference** from data (sparse and dense paths)
- [x] **Eigendecomposition** (symmetric eigendecomposition, caching)
- [x] **expQt** (matrix exponential for transition probabilities)
- [x] **Equilibrium frequencies** (pi vector, normalized)
- [x] **Exchangeability matrix** (W, symmetric, normalized by avg_transition)
- [x] **Rate scaling** (mu parameter)
- [ ] **Custom GTR from file** (`--custom-gtr` / `GTR.from_file()` in v0)
- [ ] **Random GTR generation** (`GTR.random()` in v0, used for testing)
- [ ] **GTR save to file** (`save_to_npz()` in v0)
- [ ] **Site-specific models** (multi-dimensional pi/W in v0)
- [ ] **Branch length optimization at GTR level** (v0 `optimal_t()`, `optimal_t_compressed()` - moved to partition layer in v1)
- [ ] **Sequence probability at GTR level** (v0 `prob_t()`, `prob_t_compressed()`, `prob_t_profiles()` - moved to partition layer in v1)
- [ ] **State pair compression** (v0 `state_pair()` - moved to partition layer in v1)
- [x] **JSON output** (GtrOutput struct with model type, parameters)

### 10. Alphabet System

- [x] **Nucleotide alphabet** (A, C, G, T with gap handling)
- [x] **Amino acid alphabet** (20 AAs, with/without stop codon)
- [x] **Gap-as-unknown toggle** (`treat_gap_as_unknown` per partition)
- [x] **IUPAC ambiguity codes** (nucleotide: R, Y, S, W, K, M, D, H, B, V, N, X)
- [x] **Amino acid ambiguity** (X, B, Z)
- [ ] **Alphabet guessing** from data (v0 `guess_alphabet()`: >90% nuc threshold)
- [x] **Profile maps** (ambiguity code to probability vector, in primitives layer)

### 11. Probability Distributions

- [x] **Analytical Gaussian** (PDF, product, convolution)
- [x] **Analytical Exponential** (PDF, convolution, special case a=b)
- [x] **Gaussian-Exponential convolution**
- [x] **Grid distributions** (uniform grid, linear interpolation)
- [x] **ScaledArray pattern** (normalized values + log-scale factor)
- [ ] **Unified Distribution class** (v0: wraps scipy.interpolate.interp1d)
- [ ] **Delta functions** (point masses, v0 `Distribution.delta_function()`)
- [ ] **Distribution multiplication** (v0 `Distribution.multiply()`)
- [ ] **Distribution division** (v0 `Distribution.divide()`)
- [ ] **Numerical integration** (v0 Simpson's rule, trapezoidal)
- [ ] **FFT transform** (v0 `Distribution.fft()`)
- [ ] **FWHM calculation** (full width half maximum)
- [ ] **Effective support** (threshold-based range)
- [ ] **Grid refinement** (v0 `_adjust_grid()` adaptive)
- [ ] **BranchLenInterpolator** (v0: node-specific branch length probability)
  - [ ] Input mode (Gaussian/Poisson approximation)
  - [ ] Marginal mode (from profile pairs)
  - [ ] Joint mode (from compressed state pairs)
  - [ ] Gamma rescaling (for relaxed clock)
  - [ ] Adaptive grid construction (log near zero, quadratic tails)

### 12. Representation and Partition System

- [x] **PartitionFitch** (Fitch parsimony, sparse storage)
- [x] **PartitionMarginalDense** (full probability vectors at all positions)
- [x] **PartitionMarginalSparse** (variable positions only)
- [x] **Partition traits** (PartitionMarginal, PartitionMarginalOps, PartitionCompressed, HasLogLh)
- [x] **Dense payloads** (DenseNodePartition, DenseEdgePartition, DenseSeqDis)
- [x] **Sparse payloads** (SparseNodePartition, SparseEdgePartition, MarginalSparseSeqDistribution)
- [x] **Node payloads** (NodeAncestral, NodeTimetree with time distributions)
- [x] **Edge payloads** (EdgeAncestral, EdgeTimetree with branch distributions and messages)
- [x] **Mutation composition** (edge merge: A5G + G5T = A5T with cancellation)
- [x] **Edge inversion** (swap ref/qry, reverse messages for rerooting)
- [x] **Reroot operations on partitions** (split/merge edges, update messages)
- [x] **Dense/sparse architecture** ([intentional change](../port-intentional-changes/sequence-representation-dense-sparse.md))
- [x] **Partition system** ([intentional change](../port-intentional-changes/partition-system-architecture.md))
- [x] **Graph-based tree structure** ([intentional change](../port-intentional-changes/tree-structure-graph-based.md))

### 13. Sequence Primitives

- [x] **AsciiChar** (validated u8, type-safe character)
- [x] **Seq** (Vec\<AsciiChar\> with string-like operations, concatenation, repetition)
- [x] **BitSet128/StateSet** (u128 bitfield for character sets, set operations)
- [x] **Sub** (substitution mutation: pos + ref + qry, invertible, parseable)
- [x] **InDel** (insertion/deletion: range + seq + direction, invertible)
- [x] **Composition** (character frequency tracking with mutation updates)
- [x] **find_char_ranges** (contiguous range detection for gaps, unknowns)
- [x] **compute_divs** (root-to-tip divergence calculation)

### 14. I/O Formats

**Input formats:**

- [x] FASTA (plain text)
- [x] Compressed FASTA (gzip, xz, bzip2, zstd)
- [x] Multiple FASTA files (concatenated)
- [x] Stdin input
- [x] Newick tree
- [x] Nexus tree
- [x] Phylip tree
- [x] CSV/TSV dates
- [x] CSV/TSV discrete states (mugration)
- [ ] VCF input (variant call format - not implemented in v1)
- [ ] Compressed VCF (.vcf.gz)
- [ ] Custom GTR model file

**Output formats:**

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

**v1-only formats:**

- [v1] PhyloXML
- [v1] Usher MAT (mutation-annotated tree)
- [v1] YAML serialization
- [v1] Compressed FASTA output
- [v1] Streaming readers/writers with automatic decompression

### 15. Date Parsing

- [x] Float dates (e.g. 2012.15)
- [x] ISO format (YYYY-MM-DD)
- [ ] Imprecise dates (YYYY-XX-XX)
- [x] Date ranges ([2013.2:2013.7])
- [x] Custom name/date column selection
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

## Cross-References

- [Algorithm Inventory](../port-algo-inventory/_index.md) - implementation details
- [Intentional Changes](../port-intentional-changes/_index.md) - deliberate deviations
- [Known Issues](../port-known-issues/_index.md) - bugs and missing features
- [Test Inventory](../port-test-inventory/_index.md) - test coverage
