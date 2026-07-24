# Unimplemented Algorithms

[Back to index](README.md)

Algorithms present in v0 Python that remain unported or were intentionally removed from v1 Rust. Struck-through entries marked **Ported** record capabilities that have since been implemented and should move to the relevant algorithm document when this inventory is reorganized.

---

## Joint ML

Joint maximum likelihood reconstruction finds the single most likely assignment of ancestral states across all nodes simultaneously, rather than marginalizing over alternatives at each node independently (Pupko, Pe'er, Shamir & Graur 2000). Uses traceback pointers (argmax) instead of marginalization (sum), analogous to the Viterbi algorithm for HMMs vs the forward-backward algorithm.

v0: `_ml_anc_joint()` (`#_ml_anc_joint`) in [`packages/legacy/treetime/treetime/treeanc.py#L934-L1084`](../../packages/legacy/treetime/treetime/treeanc.py#L934-L1084).
v1: retained as `MethodAncestral::Joint` only so the CLI can return an explicit removal error at [`packages/treetime/src/ancestral/pipeline.rs#L225-L232`](../../packages/treetime/src/ancestral/pipeline.rs#L225-L232). See the [intentional change](../decisions/ancestral-joint-reconstruction-removed.md).

### Algorithm

Backward pass (lines 960-1014):

- `joint_Lx[L, n_states]`: likelihood array
- `joint_Cx[L, n_states]`: argmax traceback indices (uint16)
- For each parent state j: `msg_to_parent = log(P(t))[:,j] + msg_from_children`
- `joint_Cx[:,j] = argmax(msg_to_parent)`, `joint_Lx[:,j] = max(msg_to_parent)`

Root (lines 1016-1039):

- `root.joint_Lx = msg_from_children + log(Pi)`
- Sequence chosen from normalized profile

Forward pass (lines 1043-1063):

- `node.seq_idx = choose(parent.seq_idx, node.joint_Cx.T)` - traceback via argmax pointers

Reference: <a id="cite-1"></a>[Pupko et al. 2000](https://doi.org/10.1093/oxfordjournals.molbev.a026369) [[1](#ref-1)]

---

## Local Outlier Filter

A TreeTime-specific outlier detection method that uses local parent-child timing relationships rather than global IQD-based detection. Identifies outliers by comparing each leaf's sampling date against the time implied by its local tree neighborhood, providing more sensitive detection than the global IQD method when only a few branches are anomalous.

v0: `local_filter()` (`#local_filter`), `calculate_node_timings()` (`#calculate_node_timings`), `collect_node_info()` (`#collect_node_info`), `flag_outliers()` (`#flag_outliers`) in [`packages/legacy/treetime/treetime/clock_filter_methods.py#L43-L219`](../../packages/legacy/treetime/treetime/clock_filter_methods.py#L43-L219).
v1: not ported. v1 only has the global IQD-based filter.

### Algorithm

1. `collect_node_info()` (lines 159-219): For each leaf, computes expected time from tree structure, actual sampling date, and z-score of deviation.

2. `calculate_node_timings()` (lines 108-156): Estimates internal node times using a Bayesian-like scheme: prior from parent (`t_parent + branch_length / clock_rate`), likelihood from children (weighted average of child-implied times), posterior as precision-weighted combination.

3. `flag_outliers()` (lines 71-105): Categorizes outliers by type: `date_too_early` (sampling date before expected time), `excess_mutations` (too many mutations for the time interval), `date_too_late` (sampling date after expected time).

---

## Full Covariance Matrix

Computes the N x N tip covariance matrix for phylogenetic GLS regression, where `M[i,j]` equals the sum of shared branch variances for tips i and j. v1's `ClockSet` achieves equivalent regression results via sufficient statistics propagation in O(N) without materializing the O(N^2) matrix. The full matrix computation exists in v0 for the `--covariation` mode and for producing `valid_confidence=True` (required by rate susceptibility analysis).

v0: `Cov()` (`#Cov`), `CovInv()` (`#CovInv`), `recurse()` (`#recurse`) in [`packages/legacy/treetime/treetime/treeregression.py#L125-L155`](../../packages/legacy/treetime/treetime/treeregression.py#L125-L155).
v1: not needed for basic regression (sufficient statistics approach), but missing for covariation-aware confidence estimation.

### Algorithm

`Cov()`: Recursive tree traversal accumulating shared ancestry to build the N x N matrix.

`CovInv()`: Computes the inverse recursively from child inverse blocks and rank-one corrections. This avoids calling a generic dense matrix inversion, but `full_matrix=True` still materializes an $N\times N$ result and allocates a dense matrix for every internal clade. Its cost is therefore at least $\Omega(N^2)$ and depends on topology; the repeated allocations can sum to $O(N^3)$ on a caterpillar tree. The separate sufficient-statistics regression path is the $O(N)$ algorithm.

---

## Numerical Hessian for Root Position

Augments the 2x2 Hessian (rate, intercept) with a third parameter (root position) via finite differences, enabling root position uncertainty quantification.

v0: within `find_best_root()` (`#find_best_root`) in [`packages/legacy/treetime/treetime/treeregression.py#L355-L382`](../../packages/legacy/treetime/treetime/treeregression.py#L355-L382).
v1: not ported. v1 Hessian is 2x2 (rate, intercept) only.

### Algorithm

After finding optimal root position `x*`:

1. Start with 2x2 Hessian H from rate/intercept optimization
2. Augment to 3x3 by adding root position parameter
3. Compute third row/column via finite differences (dx = +/-0.001):
   - `H[2,0] = (dCost/drate(x+dx) - dCost/drate(x-dx)) / (2*dx)`
   - `H[2,1] = (dCost/dintercept(x+dx) - dCost/dintercept(x-dx)) / (2*dx)`
   - `H[2,2] = (Cost(x+dx) - 2*Cost(x) + Cost(x-dx)) / dx^2`
4. Invert 3x3 Hessian to get covariance matrix including root position uncertainty

---

## Per-Site Rate Variation

The design document specifies a fixed rate vector $\mu^a$, where $a$ indexes alignment sites and each rate scales a shared GTR rate matrix. [`kb/_raw/sequence_evolution.md#site-specific-models`](../_raw/sequence_evolution.md#site-specific-models) This contract differs from two related models:

- Discrete-gamma among-site rate variation integrates over latent rate categories, $L_a = \sum_k w_k L_a(r_k)$, rather than assigning one known rate to each site.
- Full site-specific GTR allows parameters such as equilibrium frequencies $\pi^a$ to vary by site, so the rate-matrix eigendecomposition can also vary by site.

v0 represents per-site rates and equilibrium frequencies in `GTR_site_specific`. Its implementation recomputes transition-matrix decompositions by site even when parameters happen to permit reuse. [`packages/legacy/treetime/treetime/gtr_site_specific.py`](../../packages/legacy/treetime/treetime/gtr_site_specific.py) Exact parity therefore concerns transition probabilities and inferred outputs; v1 may share a decomposition only if it produces the same numerical result under the approved contract.

v1 stores scalar `mu` and optional `site_rates`. Its core transition matrices, dense propagation, and sparse explicit-position propagation consume the vector. Rate construction or loading, compressed fixed-position behavior, optimizer coverage, validation, serialization, CLI wiring, and end-to-end parity remain unresolved. [`packages/treetime/src/gtr/gtr.rs#L184-L193`](../../packages/treetime/src/gtr/gtr.rs#L184-L193)

Known issue: [kb/issues/M-gtr-per-site-rate-variation.md](../issues/M-gtr-per-site-rate-variation.md).

### Algorithm

For a fixed supplied rate $\mu_a$, the matrix exponential is $\exp(Q\mu_a t)$. When $Q$ is shared, an implementation may compute $Q = V\operatorname{diag}(\lambda)V^{-1}$ once and evaluate $V\operatorname{diag}(\exp(\lambda_k\mu_a t))V^{-1}$ per site. Discrete-gamma inference instead evaluates and weights multiple category-specific likelihoods for every site.

Reference: <a id="cite-2"></a>[Yang 1994](https://doi.org/10.1007/BF00160154) [[2](#ref-2)]

---

## Site-Specific GTR (Full)

Extension of GTR where equilibrium frequencies vary per alignment site, requiring per-site eigendecomposition. Standard GTR uses `Pi[n_states]`; full site-specific GTR extends to `Pi[n_states, seq_len]`. More expensive than per-site rate variation because the eigenvectors change per site.

v0: `GTR_site_specific(GTR)` (`#GTR_site_specific`) in [`packages/legacy/treetime/treetime/gtr_site_specific.py#L1-L495`](../../packages/legacy/treetime/treetime/gtr_site_specific.py#L1-L495).
v1: `GTR.is_site_specific` (`#is_site_specific`) field exists (always false); no implementation.

### Key methods

`_make_expQt_interpolator()` (`#_make_expQt_interpolator`) (lines 331-348): Pre-computes exp(Qt) on a grid of t values, then uses linear interpolation. Avoids repeated eigendecomposition during tree traversal.

`infer()` (`#infer`) (lines 208-311): Site-specific inference accumulates per-site statistics:

- `n_ija[site, state_i, state_j]`: transition counts per site
- `T_ia[site, state_i]`: total time in state i per site

---

## Stochastic Polytomy Resolution

Mutation-conditioned stochastic topology refinement as an alternative to the greedy deterministic method. The v0 generator combines mutation-removal and branch-merger events; it is not an ordinary <a id="cite-3"></a>[Kingman 1982](https://doi.org/10.1016/0304-4149(82)90011-4) [[3](#ref-3)] sampler and can leave a multifurcation unresolved.

v0: `generate_subtree()` (`#generate_subtree`) in [`packages/legacy/treetime/treetime/treetime.py#L872-L1011`](../../packages/legacy/treetime/treetime/treetime.py#L872-L1011), dispatched by `resolve_polytomies()` (`#resolve_polytomies`).
v1: not ported - v1 has greedy deterministic approach only. Known issue: [kb/issues/N-timetree-stochastic-polytomy-unimplemented.md](../issues/N-timetree-stochastic-polytomy-unimplemented.md). CLI: `--stochastic-resolve` (v0), `--greedy-resolve` (v0 inverse). v0 prints a deprecation warning for greedy mode, intending to make stochastic the default ([packages/legacy/treetime/treetime/treetime.py#L682-L685](../../packages/legacy/treetime/treetime/treetime.py#L682-L685)).

### Background

The greedy method (`_poly()` in v0, `resolve_polytomies()` in v1) repeatedly merges the pair with the highest estimated likelihood gain (<a id="cite-4"></a>[Sagulenko et al. 2018](https://doi.org/10.1093/ve/vex042) [[4](#ref-4)], Section 2.6). After a merge, the new internal node re-enters the candidate set, so repeated attachment can produce a caterpillar-like topology; the paper does not establish a general statistical bias toward that shape. The stochastic method instead samples from v0's specialized event generator. v0 intended to make stochastic the default: "Stochastic resolution will become the default in future versions" ([packages/legacy/treetime/treetime/treetime.py#L682-L685](../../packages/legacy/treetime/treetime/treetime.py#L682-L685)).

### Algorithm (`generate_subtree()`, [packages/legacy/treetime/treetime/treetime.py#L872-L1011](../../packages/legacy/treetime/treetime/treetime.py#L872-L1011))

The function models polytomy resolution as a joint mutation-coalescence process within the time window between the polytomy node and its children:

1. Sort children by `time_before_present` (most recent first)
2. Initialize `branches_alive` with the most recent child. Remaining children wait in `branches_to_come`.
3. Compute pending mutations per branch: `round(mutation_length * L)` where L = sequence length
4. Loop until two or fewer branches remain or time reaches the parent:
   - Identify branches "ready to coalesce" (zero pending mutations)
   - Compute rates: mutation rate = `gtr.mu * L * total_mutations`, coalescent rate from `merger_model.branch_merger_rate(t)` if available or dummy `2/(tmax-t)` otherwise
   - Sample waiting time: `dt = Exp(1/total_rate)` using `self.rng.exponential`
   - If a `branches_to_come` child appears in the interval, add it and restart
   - Sample event type proportional to rates (`self.rng.random`)
   - Mutation event: pick branch proportional to mutation count, decrement by one
   - Coalescent event: pick two ready branches uniformly (`self.rng.choice`), create new internal node at current time, reparent them, build `BranchLenInterpolator` for the new node
5. Remaining branches become direct children of the parent

### RNG

Uses `self.rng` (`numpy.random.default_rng(seed=rng_seed)`) from `TreeAnc.__init__()` ([packages/legacy/treetime/treetime/treeanc.py#L163](../../packages/legacy/treetime/treetime/treeanc.py#L163)). CLI flag is `--rng-seed` (v0), `--seed` (v1). Without a seed, results are non-deterministic.

### Dispatch

`resolve_polytomies(stochastic_resolve=True)` ([packages/legacy/treetime/treetime/treetime.py#L694](../../packages/legacy/treetime/treetime/treetime.py#L694)) calls `generate_subtree(n)` for each polytomy. `resolve_polytomies(stochastic_resolve=False)` calls `_poly(n)` (greedy).

---

## FFT Convolution with Delta Approximation

v0's FFT convolution with two domain-specific optimizations: delta function approximation for narrow distributions (skips FFT entirely) and linear tail extrapolation beyond the FFT valid region.

v0: `NodeInterpolator.convolve_fft()` (`#NodeInterpolator`, `#convolve_fft`) in [`packages/legacy/treetime/treetime/node_interpolator.py#L161-L267`](../../packages/legacy/treetime/treetime/node_interpolator.py#L161-L267).
v1: has FFT in treetime-ops but not the delta approximation or tail extrapolation.

### Algorithm

1. Determine grid spacing from `min(FWHM_node, FWHM_branch) / FFT_FWHM_GRID_SIZE`
2. Delta approximation (lines 185-200): If node distribution is much narrower than branch distribution (ratio < 1/fft_grid_size), skip FFT and shift branch distribution by node peak location.
3. FFT path (lines 202-240): Evaluate both distributions in plain probability space, zero-pad to `2 * raw_len` to prevent circular convolution, compute `IFFT(FFT(branch) * FFT(node))`, convert to neg-log: `res = -ln(fft_result) + peak_branch + peak_node - ln(dt)`.
4. Tail extrapolation (lines 242-260): Linearly extrapolate distribution tails beyond FFT valid region to maintain proper asymptotic behavior.

---

## Composite Simpson Convolution on an Adaptive Output Grid

Non-uniform-grid convolution that evaluates each output point by composite Simpson quadrature, then adaptively refines the output grid where interpolation error is large. It is more flexible than the uniform-grid FFT path for irregular distributions, at the cost of evaluating many convolution integrals.

v0: `NodeInterpolator.convolve()` (`#NodeInterpolator`, `#convolve`), `_evaluate_convolution()` (`#_evaluate_convolution`) in [`packages/legacy/treetime/treetime/node_interpolator.py#L268-L409`](../../packages/legacy/treetime/treetime/node_interpolator.py#L268-L409).
v1: not ported. v1 uses grid-based Riemann sum for all convolutions.

### Algorithm

1. Construct an initial non-uniform output grid: linear spacing near the predicted convolution peak and quadratic spacing in the tails.
2. At each output time $t$, evaluate $\int f(\tau)g(t-\tau)\,d\tau$ using `Distribution.integrate_simpson()`. That helper uses a fixed odd number of points on intervals split around the integrand peak; it is composite Simpson quadrature, without recursive local-error subdivision.
3. Interpolate the resulting output values, compare interpolated and directly evaluated values at candidate points, and add output-grid points where the interpolation error exceeds the tolerance. Repeat until no candidate requires refinement.

---

## FWHM Computation

Full Width at Half Maximum for discretized probability distributions. In neg-log probability space, FWHM corresponds to the region where `-log(P) < -log(P_max) + ln(2)`.

v0: `Distribution.calc_fwhm()` (`#Distribution`, `#calc_fwhm`) in [`packages/legacy/treetime/treetime/distribution.py#L24-L67`](../../packages/legacy/treetime/treetime/distribution.py#L24-L67).
v1: not implemented.

### Algorithm

1. Find peak location `x_peak` where neg-log probability is minimum
2. Find left boundary: largest x < x_peak where `-log(P(x)) = -log(P_max) + ln(2)`
3. Find right boundary: smallest x > x_peak where `-log(P(x)) = -log(P_max) + ln(2)`
4. FWHM = right_boundary - left_boundary

Used to determine appropriate grid resolution for FFT convolution and adaptive grid refinement.

---

## Branch Length Interpolator (Input Mode)

Constructs per-branch likelihood distributions over evolutionary distance from an input branch length, using a Poisson approximation for short branches and a saturation-aware Gaussian approximation for longer branches. This is v0's `branch_length_mode=input` path, which avoids sequence-based marginal likelihoods.

v0: `BranchLenInterpolator.__init__()` (`#BranchLenInterpolator`) input mode in [`packages/legacy/treetime/treetime/branch_len_interpolator.py#L64-L102`](../../packages/legacy/treetime/treetime/branch_len_interpolator.py#L64-L102).
v1: input mode currently converts each input branch length to a point distribution in [`packages/treetime/src/timetree/inference/runner.rs#L168-L192`](../../packages/treetime/src/timetree/inference/runner.rs#L168-L192); it does not implement either v0 approximation.

### Algorithm

**Short branches** (Poisson regime, lines 78-86): for candidate branch length $x$, observed mutation length $\ell$, and sequence length $L$, the negative log-likelihood up to an additive constant is $L[x-\ell\ln(x+\varepsilon)]$. The implementation uses $\varepsilon=0.01/L$ to keep the logarithm finite.

**Long branches** (Gaussian regime, lines 87-102): let $p_0=1-\sum_i\pi_i^2$, $\ell'=\ell+1/L$, and $z=\exp(\ell'/p_0)$. The branch-length variance is $\sigma^2=p_0(z-1)[z-p_0(z-1)]/L$. The implementation uses a quadratic cost near $\ell$ and a linear tail beyond a fixed cost cap to limit extreme penalties.

Transition: switches from Poisson to Gaussian at `mutation_length > 0.05`.

---

## Random GTR Generation

Testing utility that generates random GTR models from prior distributions. Useful for verifying GTR algorithms work for arbitrary (not just named) models.

v0: `GTR.random()` (`#GTR`, `#random`) in [`packages/legacy/treetime/treetime/gtr.py#L464-L489`](../../packages/legacy/treetime/treetime/gtr.py#L464-L489).
v1: not ported.

### Algorithm

1. Draw one independent `Uniform(0, 1)` value for each equilibrium-frequency entry.
2. Draw one independent `Uniform(0, 1)` value for every exchangeability-matrix entry.
3. Pass the raw arrays and caller-provided $\mu$ to `assign_rates()`, which normalizes $\pi$, clears the diagonal of $W$, rescales $W$ to TreeTime's average-rate convention, adjusts $\mu$ by the same scale, and computes the eigendecomposition. Although `assign_rates()` briefly assigns $\tfrac12(W+W^\mathsf{T})`, the next assignment overwrites that value with the scaled original matrix; random output is therefore not guaranteed to be reversible.

---

## File-Based GTR Loading

Parses a text file containing equilibrium frequencies (single-letter lines) and exchangeability parameters (letter-pair lines) to construct a custom GTR model.

v0: `GTR.from_file()` (`#GTR`, `#from_file`) in [`packages/legacy/treetime/treetime/gtr.py#L175-L233`](../../packages/legacy/treetime/treetime/gtr.py#L175-L233).
v1: not ported.

### File format

```
# Comment lines starting with #
A 0.25
C 0.25
G 0.25
T 0.25
AC 1.0
AG 2.0
AT 1.0
CG 1.0
CT 2.0
GT 1.0
```

---

## GTR Optimal Branch Length

Per-branch length optimization at the GTR level. In v1, this functionality exists in the `optimize` command via `PartitionOptimizeOps` rather than as a GTR method.

v0: `GTR.optimal_t()` (`#GTR`, `#optimal_t`), `GTR.optimal_t_compressed()` (`#optimal_t_compressed`) in [`packages/legacy/treetime/treetime/gtr.py#L788-L921`](../../packages/legacy/treetime/treetime/gtr.py#L788-L921).
v1: functionality exists in the optimize command but not as a GTR method.

### Algorithm

`optimal_t(parent_seq, child_seq)`: compress the two sequences into state pairs and multiplicities, then delegate to `optimal_t_compressed()`.

`optimal_t_compressed(state_pairs, multiplicities)`: minimize negative log-likelihood with SciPy's Brent method in $\sqrt{t}$ space, bracketed by the Hamming distance and the configured maximum branch length. An older-SciPy fallback uses bounded scalar minimization; a failed Brent result falls back to the Hamming distance.

---

## Homoplasy Scanner

Analysis pipeline that identifies recurrent mutations - sites where the same mutation occurred independently on multiple branches, indicating convergent evolution, recombination, or sequencing artifacts.

v0: `scan_homoplasies()` (`#scan_homoplasies`) in [`packages/legacy/treetime/treetime/wrappers.py#L82-L139`](../../packages/legacy/treetime/treetime/wrappers.py#L82-L139).
v1: returns an explicit not-implemented error at [`packages/treetime/src/commands/homoplasy/run.rs#L6-L8`](../../packages/treetime/src/commands/homoplasy/run.rs#L6-L8).

### Algorithm

1. Run joint ancestral reconstruction.
2. Collect non-gap, non-ambiguous substitutions by exact change $(a, p, d)$, by position, and separately on terminal branches.
3. Tabulate how often each exact substitution and each position is hit, and compare the position-hit histogram with a Poisson distribution having the same mean.
4. Rank recurrent exact substitutions by multiplicity and report the top requested entries, optionally annotated with drug-resistance metadata.

---

## N-Branches Posterior

Enhancement to the coalescent model that uses posterior probability distributions of divergence times instead of point estimates when calculating the coalescent merger rate. Produces more accurate coalescent contributions by accounting for time uncertainty.

v0: `--n-branches-posterior` flag uses posterior distributions for branch counting in merger rate computation.
v1: CLI arg declared (hidden) at [`packages/treetime/src/commands/timetree/args.rs#L184`](../../packages/treetime/src/commands/timetree/args.rs#L184), returns error at [`packages/treetime/src/commands/timetree/run.rs#L115`](../../packages/treetime/src/commands/timetree/run.rs#L115).

---

## Tree Inference from Alignment

When no input tree is provided, infer a phylogenetic tree from the alignment. v0 delegates to external tools (IQ-TREE, FastTree) or uses Bio.Phylo. Both the `timetree` and `clock` commands require a tree; these stubs guard the code path where the user omits `--tree`.

v1: `todo!()` at [`packages/treetime/src/commands/timetree/initialization.rs#L34`](../../packages/treetime/src/commands/timetree/initialization.rs#L34) and `unimplemented!()` at [`packages/treetime/src/commands/clock/run.rs#L60`](../../packages/treetime/src/commands/clock/run.rs#L60).

---

## Timetree Output: Node Dates TSV

Writes inferred node dates to TSV for downstream analysis. The output includes all tree nodes (leaves and internal).

v0: [`packages/legacy/treetime/treetime/CLI_io.py#L136-L158`](../../packages/legacy/treetime/treetime/CLI_io.py#L136-L158).
v1: `todo!()` stub at [`packages/treetime/src/commands/timetree/output/dates.rs#L20`](../../packages/treetime/src/commands/timetree/output/dates.rs#L20).
Known issue: [write_node_dates() is a todo!() stub](../issues/N-timetree-node-dates-output-unimplemented.md).

v0 format with confidence (`#node\tdate\tnumeric date\tlower bound\tupper bound`): bounds come from `get_max_posterior_region(n, fraction=0.9)`. Bad branches get `--` placeholders.

v0 format without confidence (`#node\tdate\tnumeric date`): no bounds columns.

---

## Timetree Output: Plots

Diagnostic visualizations: root-to-tip regression scatter plot and time-scaled phylogenetic tree. CLI flags `--plot-rtt` and `--plot-tree` are hidden pending implementation.

v1: `todo!()` at [`packages/treetime/src/commands/timetree/output/plots.rs#L11`](../../packages/treetime/src/commands/timetree/output/plots.rs#L11) (root-to-tip) and [`packages/treetime/src/commands/timetree/output/plots.rs#L20`](../../packages/treetime/src/commands/timetree/output/plots.rs#L20) (time tree).

---

## ~~Branch Distribution Builder~~ (Ported)

Build per-edge time likelihood distributions from partition likelihood contributions and branch-length grids, then convert their coordinates to time using the effective clock rate.

v1: `compute_branch_length_distribution()` is implemented at [`packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L23-L58`](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L23-L58) and called from [`packages/treetime/src/timetree/inference/runner.rs#L118-L127`](../../packages/treetime/src/timetree/inference/runner.rs#L118-L127). Input-branch mode remains a separate gap described above.

---

## ~~Iterative GTR for Discrete Traits~~ (Ported)

Iterative parameter estimation for discrete-trait GTR models following the expectation-maximization framework (<a id="cite-5"></a>[Dempster, Laird, and Rubin 1977](https://doi.org/10.1111/j.2517-6161.1977.tb01600.x) [[5](#ref-5)]). The E-step computes posterior joint parent-child state distributions through Felsenstein's pruning algorithm (<a id="cite-6"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[6](#ref-6)]), counting expected transitions and dwell times across all edges. The M-step re-estimates the symmetric exchangeability matrix $W$, equilibrium frequencies $\pi$, and scalar rate $\mu$ from these sufficient statistics. Rate optimization uses Brent's method with bracket validation.

v0: `reconstruct_discrete_traits()` in [`packages/legacy/treetime/treetime/wrappers.py#L785-L809`](../../packages/legacy/treetime/treetime/wrappers.py#L785-L809), `TreeAnc.infer_gtr()` in [`packages/legacy/treetime/treetime/treeanc.py#L1500-L1632`](../../packages/legacy/treetime/treetime/treeanc.py#L1500-L1632).
v1: `refine_gtr_iterative()` in [`packages/treetime/src/gtr/refinement.rs`](../../packages/treetime/src/gtr/refinement.rs). Remaining parity gap tracked in [Mugration golden master parity with v0](../issues/M-mugration-iterative-gtr.md). Full forward-backward per iteration proposed in [mugration-full-reconstruction-per-iteration](../proposals/mugration-full-reconstruction-per-iteration.md).

---

## References

- <a id="ref-1"></a>Pupko, Tal, Itsik Pe'er, Ron Shamir, and Dan Graur. 2000. "A Fast Algorithm for Joint Reconstruction of Ancestral Amino Acid Sequences." _Molecular Biology and Evolution_ 17(6):890-896. https://doi.org/10.1093/oxfordjournals.molbev.a026369 [↩](#cite-1)
- <a id="ref-2"></a>Yang, Ziheng. 1994. "Maximum Likelihood Phylogenetic Estimation from DNA Sequences with Variable Rates over Sites: Approximate Methods." _Journal of Molecular Evolution_ 39(3):306-314. https://doi.org/10.1007/BF00160154 [↩](#cite-2)
- <a id="ref-3"></a>Kingman, J. F. C. 1982. "The Coalescent." _Stochastic Processes and their Applications_ 13(3):235-248. https://doi.org/10.1016/0304-4149(82)90011-4 [↩](#cite-3)
- <a id="ref-4"></a>Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4(1):vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-4)
- <a id="ref-5"></a>Dempster, Arthur P., Nan M. Laird, and Donald B. Rubin. 1977. "Maximum Likelihood from Incomplete Data via the EM Algorithm." _Journal of the Royal Statistical Society: Series B_ 39(1):1-38. https://doi.org/10.1111/j.2517-6161.1977.tb01600.x [↩](#cite-5)
- <a id="ref-6"></a>Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _Journal of Molecular Evolution_ 17(6):368-376. https://doi.org/10.1007/BF01734359 [↩](#cite-6)
