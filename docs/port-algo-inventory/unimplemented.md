# Unimplemented Algorithms

[Back to index](_index.md)

Algorithms present in v0 Python that have not been ported to v1 Rust. Full detail for implementation planning.

---

## Joint ML

| Property    | Value                                                                                                                                                                                                              |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Type        | Well-known (Viterbi-like)                                                                                                                                                                                          |
| v0 Location | [`packages/legacy/treetime/treetime/treeanc.py#L934-L1084`](../../packages/legacy/treetime/treetime/treeanc.py#L934-L1084)                                                                                         |
| Functions   | `_ml_anc_joint()` (`#_ml_anc_joint`)                                                                                                                                                                               |
| v1 Status   | Declared in `MethodAncestral::Joint` (`#MethodAncestral`, `#Joint`) but `unimplemented!()` at [`packages/treetime/src/commands/ancestral/run.rs#L194`](../../packages/treetime/src/commands/ancestral/run.rs#L194) |
| Reference   | Pupko et al. (2000). "A fast algorithm for joint reconstruction." Mol Biol Evol, 17(6):890-896                                                                                                                     |

**Algorithm**:

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

**Difference from marginal**: Uses argmax (traceback pointers) instead of sum (marginalization). Produces single most likely reconstruction rather than posterior distribution.

---

## Local Outlier Filter

| Property    | Value                                                                                                                                                                          |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Type        | Custom (TreeTime-specific)                                                                                                                                                     |
| v0 Location | [`packages/legacy/treetime/treetime/clock_filter_methods.py#L43-L219`](../../packages/legacy/treetime/treetime/clock_filter_methods.py#L43-L219)                               |
| Functions   | `local_filter()` (`#local_filter`), `calculate_node_timings()` (`#calculate_node_timings`), `collect_node_info()` (`#collect_node_info`), `flag_outliers()` (`#flag_outliers`) |
| v1 Status   | Not ported                                                                                                                                                                     |

**Algorithm**:

1. `collect_node_info()` (lines 159-219): For each leaf, computes:
   - Expected time from tree structure
   - Actual sampling date
   - Z-score of deviation

2. `calculate_node_timings()` (lines 108-156): Estimates internal node times using a Bayesian-like scheme:
   - Prior from parent: `t_parent + branch_length / clock_rate`
   - Likelihood from children: weighted average of child-implied times
   - Posterior: precision-weighted combination

3. `flag_outliers()` (lines 71-105): Categorizes outliers:
   - `date_too_early`: sampling date before expected time
   - `excess_mutations`: too many mutations for the time interval
   - `date_too_late`: sampling date after expected time

**v0 vs v1**: v1 only has global IQD-based filter; v0 also uses local parent-child timing relationships for more sensitive outlier detection.

---

## Full Covariance Matrix

| Property    | Value                                                                                                                                  |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Standard linear algebra                                                                                                                |
| v0 Location | [`packages/legacy/treetime/treetime/treeregression.py#L125-L155`](../../packages/legacy/treetime/treetime/treeregression.py#L125-L155) |
| Functions   | `Cov()` (`#Cov`), `CovInv()` (`#CovInv`), `recurse()` (`#recurse`)                                                                     |
| v1 Status   | Not needed - sufficient statistics approach avoids materializing full matrix                                                           |

**Algorithm**:

`Cov()`: Computes N x N tip covariance matrix where `M[i,j]` = sum of shared branch variances for tips i and j. Uses recursive tree traversal to accumulate shared ancestry.

`CovInv()`: Computes inverse via Schur complement recursion:

- For a 2x2 block `[[A, B], [C, D]]`, inverse is computed from `A^{-1}` and Schur complement `D - C*A^{-1}*B`
- Recursively applied down the tree

**Note**: v1's `ClockSet` achieves equivalent results via sufficient statistics propagation in O(N) without materializing the O(N^2) covariance matrix.

---

## Rate Susceptibility Analysis

| Property    | Value                                                                                                                                                     |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Sensitivity analysis                                                                                                                                      |
| v0 Location | [`packages/legacy/treetime/treetime/clock_tree.py#L1010-L1066`](../../packages/legacy/treetime/treetime/clock_tree.py#L1010-L1066)                        |
| Functions   | `calc_rate_susceptibility()` (`#calc_rate_susceptibility`), `date_uncertainty_due_to_rate()` (`#date_uncertainty_due_to_rate`)                            |
| v1 Status   | `todo!()` at [`packages/treetime/src/commands/timetree/output/confidence.rs#L39`](../../packages/treetime/src/commands/timetree/output/confidence.rs#L39) |
| Known Issue | [--vary-rate panics with todo!()](../port-known-issues/H-timetree-vary-rate-unimplemented.md)                                                             |
| Reference   | Sagulenko et al. (2018). "TreeTime." Virus Evolution, 4(1):vex042, Section 2.5                                                                            |

**Background**: Node date uncertainty has two independent sources (Sagulenko et al. 2018). The marginal posterior from the backward/forward pass captures mutation stochasticity (source 1, implemented in v1). Clock rate uncertainty (source 2) propagates to all node dates because times scale inversely with rate. Nodes near the root have the highest sensitivity: a 10% rate error shifts the root date by 10% of the tree depth, while recent tips barely move. v0 combines both sources via quadrature sum, treating them as independent Gaussian-like contributions.

**Algorithm** (`calc_rate_susceptibility()`, [packages/legacy/treetime/treetime/clock_tree.py#L1010-L1066](../../packages/legacy/treetime/treetime/clock_tree.py#L1010-L1066)):

1. Extract `rate_std = sqrt(clock_model.cov[0,0])`. Requires either `--clock-std-dev` from the user or `valid_confidence=True` from covariation-aware regression (`--covariation` flag). Rejects negative clock rate.
2. Compute `upper = rate + rate_std`, `lower = max(0.1 * rate, rate - rate_std)`
3. Scale all edge gamma values: `gamma = gamma_orig * upper_rate / current_rate`
4. Run `make_time_tree()` with upper rate, store per-node `(upper_rate, numdate)`
5. Scale gammas for lower rate, run again, store `(lower_rate, numdate)`
6. Restore original gammas, run at central rate, store `(current_rate, numdate)`
7. Sort per-node triples by date: `[lower_date, central_date, upper_date]`

**Algorithm** (`date_uncertainty_due_to_rate()`, [packages/legacy/treetime/treetime/clock_tree.py#L1068-L1088](../../packages/legacy/treetime/treetime/clock_tree.py#L1068-L1088)):

Converts per-node `[lower, central, upper]` date triples into a Gaussian CI using `scipy.special.erfinv`:

```
nsig = [sqrt(2) * erfinv(-1 + 2*p) for p in interval]
rate_ci = [central + nsig[i] * |boundary[i] - central|]
```

**Algorithm** (`combine_confidence()`, [packages/legacy/treetime/treetime/clock_tree.py#L1090-L1101](../../packages/legacy/treetime/treetime/clock_tree.py#L1090-L1101)):

Quadrature sum of independent CI contributions:

```
lower = center - sqrt((c1_lower - center)^2 + (c2_lower - center)^2)
upper = center + sqrt((c1_upper - center)^2 + (c2_upper - center)^2)
```

Clipped to physical limits `[min_date, max_date]`.

**v0 gating**: `--confidence` alone is not sufficient. The full confidence pipeline requires EITHER `--confidence --covariation` (regression with phylogenetic covariance, auto-derives rate_std) OR `--confidence --clock-std-dev <value>` (user-specified). Without either, v0 prints a warning and disables confidence. See [--confidence flag ignored](../port-known-issues/M-timetree-confidence-flag-ignored.md).

**v0 CI method**: v0's `dates.tsv` uses `get_max_posterior_region()` (`#get_max_posterior_region`) ([packages/legacy/treetime/treetime/clock_tree.py#L1146-L1230](../../packages/legacy/treetime/treetime/clock_tree.py#L1146-L1230)), which finds the narrowest 90% highest posterior density (HPD) region around the posterior peak. v1's existing `extract_confidence_intervals()` uses `quantile(0.025, 0.975)` for an equal-tailed 95% CI. HPD regions are narrower for skewed distributions. The porting decision should clarify which approach to use.

**v1 existing infrastructure**: `combine_confidence()` (ported, unused), `Distribution::quantile()` (implemented), `ClockModel::cov()` (implemented), `EdgeTimetree.gamma` (implemented). Missing: `compute_rate_susceptibility()` body, `date_uncertainty_due_to_rate()`, per-node storage for rate-variation dates (`NodeTimetree` has no `rate_variation_dates` field).

---

## Numerical Hessian for Root Position

| Property    | Value                                                                                                                                  |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Numerical differentiation                                                                                                              |
| v0 Location | [`packages/legacy/treetime/treetime/treeregression.py#L355-L382`](../../packages/legacy/treetime/treetime/treeregression.py#L355-L382) |
| Functions   | Within `find_best_root()` (`#find_best_root`)                                                                                          |
| v1 Status   | Not ported - v1 Hessian is 2x2 (rate, intercept) only                                                                                  |

**Algorithm**:

After finding optimal root position `x*`:

1. Start with 2x2 Hessian H from rate/intercept optimization
2. Augment to 3x3 by adding root position parameter
3. Compute third row/column via finite differences (dx = +/-0.001):
   - `H[2,0] = (dCost/drate(x+dx) - dCost/drate(x-dx)) / (2*dx)`
   - `H[2,1] = (dCost/dintercept(x+dx) - dCost/dintercept(x-dx)) / (2*dx)`
   - `H[2,2] = (Cost(x+dx) - 2*Cost(x) + Cost(x-dx)) / dx^2`
4. Invert 3x3 Hessian to get covariance matrix including root position uncertainty

---

## Site-Specific GTR

| Property    | Value                                                                                                                                    |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Extension of GTR                                                                                                                         |
| v0 Location | [`packages/legacy/treetime/treetime/gtr_site_specific.py#L1-L495`](../../packages/legacy/treetime/treetime/gtr_site_specific.py#L1-L495) |
| Class       | `GTR_site_specific(GTR)` (`#GTR_site_specific`)                                                                                          |
| v1 Status   | `GTR.is_site_specific` (`#is_site_specific`) field exists (always false); no implementation                                              |

**Algorithm**:

Standard GTR has:

- `Pi[n_states]`: equilibrium frequencies (single vector)
- `W[n_states, n_states]`: exchangeability matrix
- `mu`: scalar rate

Site-specific GTR extends to:

- `Pi[n_states, seq_len]`: per-site equilibrium frequencies
- `mu[seq_len]`: per-site rates

Key methods:

`_make_expQt_interpolator()` (`#_make_expQt_interpolator`) (lines 331-348): Pre-computes exp(Qt) on a grid of t values, then uses linear interpolation. Avoids repeated eigendecomposition during tree traversal.

`infer()` (`#infer`) (lines 208-311): Site-specific inference accumulates:

- `n_ija[site, state_i, state_j]`: transition counts per site
- `T_ia[site, state_i]`: total time in state i per site

---

## Stochastic Polytomy Resolution

| Property    | Value                                                                                                                        |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------- |
| Type        | Monte Carlo / coalescent simulation                                                                                          |
| v0 Location | [`packages/legacy/treetime/treetime/treetime.py#L872-L1011`](../../packages/legacy/treetime/treetime/treetime.py#L872-L1011) |
| Functions   | `generate_subtree()` (`#generate_subtree`), `resolve_polytomies()` (`#resolve_polytomies`)                                   |
| v1 Status   | Not ported - v1 only has greedy deterministic approach                                                                       |
| Known Issue | [Stochastic polytomy resolution not implemented](../port-known-issues/N-timetree-stochastic-polytomy-unimplemented.md)       |
| CLI         | `--stochastic-resolve` (v0), `--greedy-resolve` (v0 inverse). v0 prints deprecation warning for greedy mode.                 |
| Reference   | Kingman (1982). "The coalescent." Stochastic Processes and Applications, 13(3):235-248                                       |

**Background**: The greedy method (`_poly()` in v0, `resolve_polytomies()` in v1) always merges the pair with the highest likelihood gain. This biases toward caterpillar-like (comb) topologies because after the first merge creates a new internal node, subsequent merges preferentially attach to it (Sagulenko et al. 2018, Section 2.6). The stochastic method samples resolutions from the Kingman coalescent process, producing tree shapes consistent with population dynamics. v0 intended to make stochastic the default: "Stochastic resolution will become the default in future versions" ([packages/legacy/treetime/treetime/treetime.py#L682-L685](../../packages/legacy/treetime/treetime/treetime.py#L682-L685)).

**Algorithm** (`generate_subtree()`, [packages/legacy/treetime/treetime/treetime.py#L872-L1011](../../packages/legacy/treetime/treetime/treetime.py#L872-L1011)):

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
   - **Mutation event**: pick branch proportional to mutation count, decrement by one
   - **Coalescent event**: pick two ready branches uniformly (`self.rng.choice`), create new internal node at current time, reparent them, build `BranchLenInterpolator` for the new node
5. Remaining branches become direct children of the parent

**RNG**: Uses `self.rng` (`numpy.random.default_rng(seed=rng_seed)`) from `TreeAnc.__init__()` ([packages/legacy/treetime/treetime/treeanc.py#L163](../../packages/legacy/treetime/treetime/treeanc.py#L163)). CLI flag is `--rng-seed` (v0), `--seed` (v1). Without a seed, results are non-deterministic.

**Dispatch**: `resolve_polytomies(stochastic_resolve=True)` ([packages/legacy/treetime/treetime/treetime.py#L694](../../packages/legacy/treetime/treetime/treetime.py#L694)) calls `generate_subtree(n)` for each polytomy. `resolve_polytomies(stochastic_resolve=False)` calls `_poly(n)` (greedy).

---

## FFT Convolution with Delta Approximation

| Property    | Value                                                                                                                                        |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known (FFT) + domain heuristic                                                                                                          |
| v0 Location | [`packages/legacy/treetime/treetime/node_interpolator.py#L161-L267`](../../packages/legacy/treetime/treetime/node_interpolator.py#L161-L267) |
| Functions   | `NodeInterpolator.convolve_fft()` (`#NodeInterpolator`, `#convolve_fft`)                                                                     |
| v1 Status   | v1 has FFT in treetime-ops but not the delta approximation or tail extrapolation                                                             |

**Algorithm**:

1. Determine grid spacing from `min(FWHM_node, FWHM_branch) / FFT_FWHM_GRID_SIZE`

2. Delta approximation (lines 185-200): If node distribution is much narrower than branch distribution (ratio < 1/fft_grid_size), treat node as delta function:
   - Skip FFT entirely
   - Just shift branch distribution by node peak location

3. FFT path (lines 202-240):
   - Evaluate both distributions in plain probability space
   - Zero-pad to `2 * raw_len` to prevent circular convolution
   - `fft_result = IFFT(FFT(branch) * FFT(node))`
   - Convert to neg-log: `res = -ln(fft_result) + peak_branch + peak_node - ln(dt)`

4. Tail extrapolation (lines 242-260): Linearly extrapolate distribution tails beyond FFT valid region to maintain proper asymptotic behavior

---

## Adaptive Simpson's Rule Convolution

| Property    | Value                                                                                                                                        |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known (adaptive quadrature)                                                                                                             |
| v0 Location | [`packages/legacy/treetime/treetime/node_interpolator.py#L268-L409`](../../packages/legacy/treetime/treetime/node_interpolator.py#L268-L409) |
| Functions   | `NodeInterpolator.convolve()` (`#NodeInterpolator`, `#convolve`), `_evaluate_convolution()` (`#_evaluate_convolution`)                       |
| v1 Status   | Not ported - v1 uses grid-based Riemann sum                                                                                                  |

**Algorithm**:

1. Construct non-uniform grid (lines 280-310):
   - Linear spacing near distribution peaks
   - Quadratic spacing in tails
   - Adaptive refinement based on local curvature

2. At each output grid point t, evaluate convolution integral (lines 320-360):
   - Construct integrand `f(tau) * g(t - tau)`
   - Apply Simpson's rule: `integral ≈ (h/3) * [f(a) + 4*f(m) + f(b)]`
   - Subdivide if error estimate exceeds tolerance

3. Refinement (lines 370-400):
   - Compare interpolated vs computed values
   - Add grid points where interpolation error is large
   - Iterate until convergence

**Purpose**: More accurate than FFT for irregular distributions; slower but handles edge cases better.

---

## FWHM Computation

| Property    | Value                                                                                                                          |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------ |
| Type        | Standard statistics                                                                                                            |
| v0 Location | [`packages/legacy/treetime/treetime/distribution.py#L24-L67`](../../packages/legacy/treetime/treetime/distribution.py#L24-L67) |
| Functions   | `Distribution.calc_fwhm()` (`#Distribution`, `#calc_fwhm`)                                                                     |
| v1 Status   | Not implemented                                                                                                                |

**Algorithm**:

In neg-log probability space, FWHM corresponds to region where `-log(P) < -log(P_max) + ln(2)`:

1. Find peak location `x_peak` where neg-log probability is minimum
2. Find left boundary: largest x < x_peak where `-log(P(x)) = -log(P_max) + ln(2)`
3. Find right boundary: smallest x > x_peak where `-log(P(x)) = -log(P_max) + ln(2)`
4. FWHM = right_boundary - left_boundary

**Purpose**: Used to determine appropriate grid resolution for FFT convolution and adaptive grid refinement.

---

## Branch Length Interpolator (Input Mode)

| Property    | Value                                                                                                                                                                                                                      |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Poisson/Gaussian approximation                                                                                                                                                                                             |
| v0 Location | [`packages/legacy/treetime/treetime/branch_len_interpolator.py#L64-L102`](../../packages/legacy/treetime/treetime/branch_len_interpolator.py#L64-L102)                                                                     |
| Functions   | `BranchLenInterpolator.__init__()` (`#BranchLenInterpolator`) input mode                                                                                                                                                   |
| v1 Status   | Partially in [`packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs`](../../packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs) but without Poisson/Gaussian switching |

**Algorithm**:

Short branches (Poisson regime, lines 78-86):

- `log_prob = -k*L + l*L*log(k)` where k = branch_length, l = mutation count, L = seq_length
- Valid when expected mutations << sequence length

Long branches (Gaussian regime, lines 87-102):

- Mean: `mu = mutation_count / clock_rate`
- Variance: `sigma^2 = mutation_count * (1 + overdispersion) / clock_rate^2`
- `log_prob = -0.5 * ((t - mu) / sigma)^2`
- Accounts for substitution saturation

Transition: Switches from Poisson to Gaussian at `mutation_length > 0.05`.

---

## Random GTR Generation

| Property    | Value                                                                                                            |
| ----------- | ---------------------------------------------------------------------------------------------------------------- |
| Type        | Testing utility                                                                                                  |
| v0 Location | [`packages/legacy/treetime/treetime/gtr.py#L464-L489`](../../packages/legacy/treetime/treetime/gtr.py#L464-L489) |
| Functions   | `GTR.random()` (`#GTR`, `#random`)                                                                               |
| v1 Status   | Not ported                                                                                                       |

**Algorithm**:

1. Sample equilibrium frequencies: `pi = Dirichlet(alpha=1)` (uniform on simplex)
2. Sample symmetric exchangeability matrix: `W[i,j] = W[j,i] = Exponential(1)`
3. Normalize: `mu = 1 / sum(pi_i * sum_j(W[i,j] * pi_j))`
4. Build rate matrix Q and compute eigendecomposition

**Purpose**: Used in unit tests to verify GTR algorithms work for arbitrary (not just named) models.

---

## File-Based GTR Loading

| Property    | Value                                                                                                            |
| ----------- | ---------------------------------------------------------------------------------------------------------------- |
| Type        | I/O utility                                                                                                      |
| v0 Location | [`packages/legacy/treetime/treetime/gtr.py#L175-L233`](../../packages/legacy/treetime/treetime/gtr.py#L175-L233) |
| Functions   | `GTR.from_file()` (`#GTR`, `#from_file`)                                                                         |
| v1 Status   | Not ported                                                                                                       |

**Algorithm**:

Parses text file with format:

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

Frequencies (single letter) and exchangeabilities (letter pairs) are parsed and used to construct GTR model.

---

## GTR Optimal Branch Length

| Property    | Value                                                                                                            |
| ----------- | ---------------------------------------------------------------------------------------------------------------- |
| Type        | Optimization                                                                                                     |
| v0 Location | [`packages/legacy/treetime/treetime/gtr.py#L788-L921`](../../packages/legacy/treetime/treetime/gtr.py#L788-L921) |
| Functions   | `GTR.optimal_t()` (`#GTR`, `#optimal_t`), `GTR.optimal_t_compressed()` (`#optimal_t_compressed`)                 |
| v1 Status   | Functionality exists in `optimize` command but not as GTR method                                                 |

**Algorithm**:

`optimal_t(parent_seq, child_seq)`:

1. Compute profile from sequences
2. Evaluate log-likelihood on grid of t values
3. Find argmax via golden section or Brent

`optimal_t_compressed(state_pairs, multiplicities)`:

- Same but uses compressed representation of aligned sequence pairs
- `state_pairs[k] = (parent_state, child_state)`
- `multiplicities[k]` = count of this pair in alignment

---

## Homoplasy Scanner

| Property    | Value                                                                                                                                  |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Analysis pipeline                                                                                                                      |
| v0 Location | [`packages/legacy/treetime/treetime/wrappers.py#L82-L139`](../../packages/legacy/treetime/treetime/wrappers.py#L82-L139)               |
| Functions   | `scan_homoplasies()` (`#scan_homoplasies`)                                                                                             |
| v1 Status   | `unimplemented!()` at [`packages/treetime/src/commands/homoplasy/run.rs#L5`](../../packages/treetime/src/commands/homoplasy/run.rs#L5) |

**Algorithm**:

1. Run ancestral reconstruction (marginal or joint)
2. Collect mutations per site across all branches
3. Identify homoplasies: sites where the same mutation occurred independently on multiple branches
4. Compute statistics: mutation count per site, consistency index, retention index
5. Output per-site mutation table with branch annotations

**Purpose**: Identifies recurrent mutations that arose independently on separate lineages, indicating convergent evolution, recombination, or sequencing artifacts.

---

## N-Branches Posterior

| Property    | Value                                                                                                                                                                                                                                                                                  |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Enhancement to coalescent model                                                                                                                                                                                                                                                        |
| v0 Location | v0 `--n-branches-posterior` flag uses posterior distributions of divergence times for estimating number of branches in coalescent merger rate                                                                                                                                          |
| v1 Status   | CLI arg declared (hidden) at [`packages/treetime/src/commands/timetree/args.rs#L184`](../../packages/treetime/src/commands/timetree/args.rs#L184), returns error at [`packages/treetime/src/commands/timetree/run.rs#L115`](../../packages/treetime/src/commands/timetree/run.rs#L115) |

**Purpose**: Uses posterior probability distributions of divergence times instead of point estimates when calculating the coalescent merger rate. Produces more accurate coalescent contributions by accounting for time uncertainty.

---

## Tree Inference from Alignment

| Property    | Value                                                                                                                                                                                                                                                                                    |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Tree building                                                                                                                                                                                                                                                                            |
| v0 Location | v0 delegates to external tools (IQ-TREE, FastTree) or uses Bio.Phylo                                                                                                                                                                                                                     |
| v1 Status   | `todo!()` at [`packages/treetime/src/commands/timetree/initialization.rs#L34`](../../packages/treetime/src/commands/timetree/initialization.rs#L34) and `unimplemented!()` at [`packages/treetime/src/commands/clock/run.rs#L60`](../../packages/treetime/src/commands/clock/run.rs#L60) |

**Purpose**: When no input tree is provided, infer a phylogenetic tree from the alignment. Both the `timetree` and `clock` commands require a tree; these stubs guard the code path where the user omits `--tree`.

---

## Timetree Output: Node Dates TSV

| Property    | Value                                                                                                                                           |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | I/O utility                                                                                                                                     |
| v0 Location | [`packages/legacy/treetime/treetime/CLI_io.py#L136-L158`](../../packages/legacy/treetime/treetime/CLI_io.py#L136-L158)                          |
| v1 Status   | `todo!()` at [`packages/treetime/src/commands/timetree/output/dates.rs#L20`](../../packages/treetime/src/commands/timetree/output/dates.rs#L20) |
| Known Issue | [write_node_dates() is a todo!() stub](../port-known-issues/N-timetree-node-dates-output-unimplemented.md)                                      |

**Purpose**: Write inferred node dates to TSV for downstream analysis.

**v0 format** (with confidence): `#node\tdate\tnumeric date\tlower bound\tupper bound`. Bounds from `get_max_posterior_region(n, fraction=0.9)`. Bad branches get `--` placeholders.

**v0 format** (without confidence): `#node\tdate\tnumeric date`.

---

## Timetree Output: Plots

| Property  | Value                                                                                                                                                                                                                                                                                                            |
| --------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type      | Visualization                                                                                                                                                                                                                                                                                                    |
| v1 Status | `todo!()` at [`packages/treetime/src/commands/timetree/output/plots.rs#L11`](../../packages/treetime/src/commands/timetree/output/plots.rs#L11) (root-to-tip) and [`packages/treetime/src/commands/timetree/output/plots.rs#L20`](../../packages/treetime/src/commands/timetree/output/plots.rs#L20) (time tree) |

**Purpose**: Generate diagnostic visualizations - root-to-tip regression scatter plot and time-scaled phylogenetic tree. CLI flags `--plot-rtt` and `--plot-tree` are hidden pending implementation.

---

## Branch Distribution Builder

| Property  | Value                                                                                                                                                                                                         |
| --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type      | Distribution construction                                                                                                                                                                                     |
| v1 Status | `todo!()` at [`packages/treetime/src/commands/timetree/inference/branch_distributions.rs#L39`](../../packages/treetime/src/commands/timetree/inference/branch_distributions.rs#L39) - function body is a stub |

**Purpose**: Build per-edge time prior distributions from branch lengths or marginal reconstruction messages. The function signature and types exist; the implementation body is a stub.

---

## Nexus File Reading

| Property  | Value                                                                                                                              |
| --------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| Type      | I/O utility                                                                                                                        |
| v1 Status | `unimplemented!()` at [`packages/treetime-cli/src/convert/convert.rs#L90`](../../packages/treetime-cli/src/convert/convert.rs#L90) |

**Purpose**: Parse Nexus tree format. The `convert` CLI command supports multiple tree formats but Nexus reading is not yet implemented.

---

## Iterative GTR for Discrete Traits

| Property    | Value                                                                                                                                                                                                                                                                                                                                                    |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Iterative parameter estimation                                                                                                                                                                                                                                                                                                                           |
| v0 Location | [`packages/legacy/treetime/treetime/wrappers.py#L785-L809`](../../packages/legacy/treetime/treetime/wrappers.py#L785-L809)                                                                                                                                                                                                                               |
| Functions   | `reconstruct_discrete_traits()` (`#reconstruct_discrete_traits`), `TreeAnc.infer_gtr()` (`#infer_gtr`) at [`treeanc.py#L1500-L1632`](../../packages/legacy/treetime/treetime/treeanc.py#L1500-L1632), `TreeAnc.optimize_gtr_rate()` (`#optimize_gtr_rate`) at [`treeanc.py#L1679-L1708`](../../packages/legacy/treetime/treetime/treeanc.py#L1679-L1708) |
| v1 Status   | Not ported. Mugration uses uniform GTR at [`packages/treetime/src/commands/mugration/run.rs#L124-L129`](../../packages/treetime/src/commands/mugration/run.rs#L124-L129)                                                                                                                                                                                 |

**Algorithm**:

After initial marginal ancestral reconstruction with `infer_gtr=True`:

1. `infer_gtr()` counts state transitions across all branches using the current marginal reconstruction. Estimates new equilibrium frequencies and symmetric exchangeability matrix from the observed transition/dwell-time statistics.
2. `optimize_gtr_rate()` optimizes the scalar rate `mu` via `scipy.optimize.minimize_scalar` (Brent), minimizing negative log-likelihood of the observed branch lengths under the current GTR model.
3. Steps 1-2 repeat for 5 iterations (default `iterations` parameter in `reconstruct_discrete_traits()`).
4. Final `infer_ancestral_sequences(infer_gtr=False)` reconstructs with the refined model.

The iterative refinement produces non-uniform equilibrium frequencies that reflect actual trait prevalence in the data. For mugration (geographic traits), this means common locations receive higher prior weight. Without iterative GTR, v1 assigns uniform prior weight to all locations, causing the argmax to differ at ambiguous internal nodes where the phylogeographic signal is weak.

**Impact**: Golden master tests show 4/6 datasets diverge from v0 at internal nodes (dengue, tb, rsv, mpox). 2/6 datasets agree despite uniform rates (zika, lassa) because their phylogeographic signal is strong enough to overwhelm the prior difference.

See [Iterative GTR inference not implemented for mugration](../port-known-issues/M-mugration-iterative-gtr.md) for the known issue entry.
