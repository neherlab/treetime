# Branch length optimization algorithms: v0, v1, and the field

This report covers the per-edge optimizer, the outer convergence loop, and the timetree integration. For the initial guess / starting point, see [optimize-initial-branch-length-estimation.md](optimize-initial-branch-length-estimation.md).

## The three-layer optimization architecture

Phylogenetic branch length optimization has three nested layers. This structure is the same across all tools (<a id="cite-12"></a>[Stamatakis, 2014](https://doi.org/10.1093/bioinformatics/btu033) [[12](#ref-12)] RAxML-NG, <a id="cite-13"></a>[Nguyen et al., 2015](https://doi.org/10.1093/molbev/msu300) [[13](#ref-13)], <a id="cite-14"></a>[Minh et al., 2020](https://doi.org/10.1093/molbev/msaa015) [[14](#ref-14)] IQ-TREE, PhyML, <a id="cite-15"></a>[Sagulenko et al., 2018](https://doi.org/10.1093/ve/vex042) [[15](#ref-15)] TreeTime) and is a form of coordinate descent (block coordinate ascent on the log-likelihood).

**Inner layer: per-edge optimizer.** Given fixed ancestral state distributions at both endpoints, find the branch length `t` that maximizes the conditional log-likelihood for one edge. This is a 1D scalar optimization problem.

**Middle layer: coordinate sweep.** Iterate over all edges in the tree, optimizing each one with the inner layer. One complete sweep over all edges is one "round" of optimization.

**Outer layer: alternating optimization.** Alternate between (a) fixing branch lengths and recomputing ancestral state distributions (marginal reconstruction), and (b) fixing distributions and optimizing all branch lengths (coordinate sweep). Repeat until convergence.

The outer layer is an EM-like algorithm: marginal reconstruction is the E-step (compute expected sufficient statistics), branch length optimization is the M-step (maximize parameters given statistics). This guarantees monotonic likelihood increase per cycle <a id="cite-8"></a>[Hobolth & Yoshida, 2005](https://arxiv.org/abs/q-bio/0511034) [[8](#ref-8)].

## Layer 1: per-edge optimization

### The likelihood function

The per-site conditional likelihood derives from <a id="cite-10"></a>[Felsenstein, 1981](https://doi.org/10.1007/BF01734359) [[10](#ref-10)] pruning. For a single edge connecting parent node `p` to child node `c`:

```
L_i(t) = sum_{a,b} msg_child_i(a) * exp(Q*t)_{ab} * msg_parent_i(b)
```

where `msg_child` is the outgroup message (evidence from the rest of the tree about the parent's state) and `msg_parent` is the subtree message (evidence from below about the child's state). The eigendecomposition `Q = V * diag(lambda) * V^{-1}` enables efficient evaluation:

```
L_i(t) = sum_c k_{ic} * exp(lambda_c * t)
```

where the coefficients `k_{ic} = (msg_child.dot(V))_c * (msg_parent.dot(V_inv.T))_c` are precomputed once per edge and reused across Newton iterations. The first and second derivatives follow by multiplying by `lambda_c` and `lambda_c^2`:

```
L'(t)  = sum_i [sum_c k_{ic} * lambda_c * exp(lambda_c * t)] / L_i(t)
L''(t) = sum_i [sum_c k_{ic} * lambda_c^2 * exp(lambda_c * t)] / L_i(t) - [L'_i(t)]^2
```

Analytical derivatives via the eigendecomposition are much faster than automatic differentiation <a id="cite-4"></a>[Fourment et al., 2022](https://arxiv.org/abs/2211.02168) [[4](#ref-4)]. Per-site likelihood updates can be further reduced to O(log n) with the LvD algorithm <a id="cite-5"></a>[Bryant et al., 2025](https://arxiv.org/abs/2601.19064) [[5](#ref-5)].

### Unimodality

<a id="cite-2"></a>[Dinh & Matsen, 2015](https://arxiv.org/abs/1507.03647) [[2](#ref-2)] prove the per-edge log-likelihood has at most one stationary point (unimodal) under JC69, F81, and all binary symmetric models. For K2P (Kimura 2-parameter) they construct examples with two local maxima. The space of rescaled 1D likelihood functions under K2P is dense in continuous functions on `[0, inf)` - any shape is possible.

Practical implication: for JC69/F81, any root-finding method converges to the global optimum. For K80, HKY85, GTR, safeguards are needed (curvature checks, fallback to bracket-based methods, grid search).

### Implementation comparison

| Aspect                         | v0                           | v1                           | RAxML-NG                     | IQ-TREE                               | PhyML                    |
| ------------------------------ | ---------------------------- | ---------------------------- | ---------------------------- | ------------------------------------- | ------------------------ |
| **Method**                     | Brent in sqrt(t) space       | Newton-Raphson               | Newton-Raphson               | NR + bisection hybrid                 | Cubic spline on f'       |
| **Derivatives**                | None (function only)         | 1st + 2nd analytical         | 1st + 2nd analytical         | 1st + 2nd analytical                  | 1st only (2 points)      |
| **Non-concave fallback**       | N/A (Brent is bracket-based) | Grid search (100 pts)        | Step = -f/\|df\|             | Bisection step                        | Bracket shrinking        |
| **Optimizer failure fallback** | Hamming distance             | Keep current value           | Revert (SAFE mode)           | Brent (`--no-opt-by-newton`)          | Brent (available)        |
| **Max inner iterations**       | scipy default                | 10                           | 30                           | 100                                   | 1020                     |
| **Branch bounds**              | [0, 4.0]                     | [0, inf)                     | [1e-6, 100]                  | [1e-6, 10]                            | [1e-8, 100]              |
| **Near-zero handling**         | sqrt(t) reparameterization   | Zero-branch derivative check | z=exp(-t) reparameterization | Clamp to min_branch                   | Geometric bracket search |
| **Regularization**             | exp(t^4/10000) penalty       | None                         | None                         | Divergence guard (revert if >95% max) | None                     |

### v0 implementation detail

`GTR.optimal_t_compressed(profiles=True)` uses <a id="cite-1"></a>[Brent, 1973](https://maths-people.anu.edu.au/~brent/pd/rpb011i.pdf) [[1](#ref-1)] scalar minimization. Code at [packages/legacy/treetime/treetime/gtr.py#L816-L920](../../packages/legacy/treetime/treetime/gtr.py#L816-L920):

```python
def _neg_prob(t, seq_pair, multiplicity):
    # Optimize in sqrt(t) space: input is s, actual branch length is s^2
    res = -1.0 * self.prob_t_profiles(seq_pair, multiplicity, t**2, return_log=True)
    return res + np.exp(t**4 / 10000)  # regularization penalty

hamming_distance = 1 - np.sum(multiplicity * np.sum(seq_pair[0] * seq_pair[1], axis=1)) / np.sum(multiplicity)

opt = minimize_scalar(_neg_prob,
    bracket=[-np.sqrt(MAX_BRANCH_LENGTH), np.sqrt(hamming_distance), np.sqrt(MAX_BRANCH_LENGTH)],
    args=(seq_pair, multiplicity), tol=tol, method='brent')

new_len = opt['x'] ** 2   # convert back from sqrt-space

if opt['success'] != True:
    new_len = hamming_distance  # fallback
```

The sqrt(t) reparameterization serves two purposes: (1) enforces `t >= 0` via `t = s^2`, and (2) stretches the likelihood surface near zero where `dL/dt` can be steep, improving Brent's convergence. The regularization `exp(s^4/10000) = exp(t^2/10000)` penalizes long branches (grows superexponentially beyond `t ~ 10`).

### v1 implementation detail

`run_optimize_mixed()` at [packages/treetime/src/commands/optimize/optimize_unified.rs#L216-L289](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L216-L289):

```rust
// Zero-branch short-circuit: check derivative sign at t=0
if is_zero_branch_optimal(&contributions) {
    edge.set_branch_length(Some(0.0));
    return Ok(());
}

// Newton-Raphson
if metrics.second_derivative < 0.0 {
    new_branch_length = branch_length - clamp(derivative / second_derivative, -1.0, branch_length);
    // iterate until |delta| <= 0.001 * bl, max 10 iterations
}

// Grid search fallback (when Hessian is non-negative)
else {
    let branch_lengths = Array1::linspace(0.1 * one_mutation, 1.5 * branch_length + one_mutation, 100);
    // select argmax likelihood
}
```

The Newton step is clamped to `[-1.0, current_bl]` for non-negativity. The zero-branch check evaluates `sum_i (sum_c k_{ic} * lambda_c) / (sum_c k_{ic})` - the derivative of log-likelihood at `t=0` as a weighted average of eigenvalues.

### v1 gaps relative to the field

1. **No hard branch length upper bound.** All other tools cap at 4-100. v1's grid search upper bound `1.5*t + 1/L` is data-dependent with no absolute cap.
2. **No reparameterization.** v0 uses sqrt(t), RAxML uses z=exp(-t). v1 operates in natural t-space with step clamping.
3. **No per-branch likelihood check.** RAxML's SAFE mode reverts a branch if likelihood decreased. v1 accepts the Newton/grid result unconditionally.
4. **Grid search range depends on initial guess.** If the initial guess underestimates (e.g., hard Hamming gives 0 at uncertain edges), the grid is too narrow.

## Layer 2: coordinate sweep

The coordinate sweep iterates over all edges, optimizing each one independently. In v0 this happens inside `optimize_tree_marginal`. In v1, `run_optimize_mixed` iterates over `graph.get_edges()`.

Both v0 and v1 do a single complete sweep per outer iteration (all edges once). Some ML tools (RAxML-NG) do multiple sweeps ("smoothings") per outer iteration - up to 32.

The coordinate sweep guarantees monotonic improvement: each per-edge optimization can only increase (or maintain) the total log-likelihood, because the likelihood decomposes as a product over edges.

## Layer 3: outer loop (alternating optimization)

### v0 outer loop

`optimize_tree_marginal()` at [packages/legacy/treetime/treetime/treeanc.py#L1297-L1360](../../packages/legacy/treetime/treetime/treeanc.py#L1297-L1360):

```python
for i in range(max_iter):
    if infer_gtr:
        self.infer_gtr(marginal=True)
        self.infer_ancestral_sequences(marginal=True)

    tol = 1e-8 + 0.01 ** (i + 1)              # progressive tolerance
    for n in self.tree.find_clades():
        new_val = self.optimal_marginal_branch_length(n, tol=tol)
        # exponential damping
        update_val = new_val * (1 - 0.75 ** (i + 1)) + n.branch_length * 0.75 ** (i + 1)
        n.branch_length = update_val

    self.infer_ancestral_sequences(marginal=True)   # recompute profiles

    if deltaLH < LHtol:   # LHtol = 0.1
        break
```

Three stabilization mechanisms:

1. **Exponential damping**: `alpha(i) = 1 - 0.75^(i+1)`. Iteration 0: 25% new, 75% old. Converges to full update. Prevents oscillation by making early steps conservative.
2. **Progressive tolerance**: Brent tolerance tightens from 0.01 to 1e-8 over iterations. Early iterations use coarse optimization, later iterations refine.
3. **Convergence on log-likelihood**: stops when `deltaLH < 0.1`.

### v1 outer loop

`run_optimize()` at [packages/treetime/src/commands/optimize/run.rs#L126-L146](../../packages/treetime/src/commands/optimize/run.rs#L126-L146):

```rust
for i in 0..*max_iter {
    let sparse_lh = update_marginal(&graph, &sparse_partitions)?;
    let dense_lh = update_marginal(&graph, &dense_partitions)?;
    let total_lh = sparse_lh + dense_lh;

    if (total_lh - lh_prev).abs() < dp.abs() {
        break;
    }

    run_optimize_mixed(&graph, &mixed_partitions)?;
    lh_prev = total_lh;
}
```

No stabilization mechanisms. Full Newton-optimal branch lengths replace previous values each iteration. Known to cause oscillation (documented in [M-optimize-oscillation-no-damping.md](../port-known-issues/M-optimize-oscillation-no-damping.md)).

### Comparison with ML tools

| Tool        | Outer damping          | Stabilization strategy                                      | Max outer iters   | Convergence         |
| ----------- | ---------------------- | ----------------------------------------------------------- | ----------------- | ------------------- |
| RAxML-NG    | None                   | Per-branch step clamping (`xmax/iters`), SAFE mode (revert) | 32                | LH change < epsilon |
| IQ-TREE     | None                   | Per-branch bisection fallback, divergence guard             | context-dependent | LH change < 0.001   |
| PhyML       | None                   | Per-branch bracket shrinking                                | until convergence | LH change < 1e-3    |
| TreeTime v0 | Exponential 0.75^(i+1) | Progressive tolerance, regularization penalty               | 10                | LH change < 0.1     |
| TreeTime v1 | None                   | None                                                        | user-specified    | LH change < dp      |

ML tools do not use outer-loop damping. They rely on per-branch safeguards to prevent individual edges from overshooting. TreeTime v0's exponential damping is unusual in the field.

### Convergence theory

<a id="cite-6"></a>[Clancy et al., 2025](https://arxiv.org/abs/2507.22038) [[6](#ref-6)] prove that coordinate maximization of the phylogenetic log-likelihood converges exponentially fast in the strong-signal regime (deep inside the Kesten-Stigum reconstruction threshold, balanced tree, poly(n) alignment columns). Outside this regime, no general convergence rate guarantee exists.

<a id="cite-7"></a>[Mossel et al., 2008](https://arxiv.org/abs/0802.0914) [[7](#ref-7)] prove that joint (max-product) ancestral reconstruction is statistically inconsistent for branch length estimation - it systematically shrinks short branches. Marginal (sum-product) reconstruction is consistent. v1 correctly uses the two-pass marginal algorithm <a id="cite-11"></a>[Pupko et al., 2000](https://doi.org/10.1093/oxfordjournals.molbev.a026369) [[11](#ref-11)].

### Advanced: Anderson acceleration

<a id="cite-9a"></a>[Henderson & Varadhan, 2018](https://arxiv.org/abs/1803.06673) [[9](#ref-9)] show that Anderson acceleration with damping and periodic restarts dramatically speeds up EM-like iterations. It is a drop-in replacement for the fixed-point iteration (no model-specific tuning needed). No phylogenetic tool currently uses Anderson acceleration. This could be a v1 innovation for the outer loop.

## Timetree integration

### v0 timetree pipeline

`TreeTime.run()` at [packages/legacy/treetime/treetime/treetime.py#L220-L350](../../packages/legacy/treetime/treetime/treetime.py#L220-L350):

```
# Pre-loop: two rounds of ML branch optimization
if branch_length_mode != 'input':
    optimize_tree(infer_gtr=True, max_iter=1)   # line 243: 1 round with GTR inference
    ...
    optimize_tree(max_iter=1)                     # line 266: 1 round after reroot

make_time_tree()                                  # first time inference

# Iteration loop:
for niter in range(max_iter):
    if polytomies_resolved:
        optimize_tree(max_iter=0)                 # line 334: reconstruction only, no optimization
    if need_new_time_tree:
        make_time_tree()
        infer_ancestral_sequences()
    else:
        infer_ancestral_sequences()
        make_time_tree()
```

Key points:

- ML branch optimization runs **before** the first `make_time_tree()`, not inside the iteration loop
- Inside the loop, `optimize_tree(max_iter=0)` only runs reconstruction (no branch optimization) and only after polytomy resolution
- The iteration loop alternates between ancestral reconstruction and time inference, not between reconstruction and ML branch optimization

### v1 timetree pipeline

`run_timetree_estimation()` at [packages/treetime/src/commands/timetree/run.rs#L35-L324](../../packages/treetime/src/commands/timetree/run.rs#L35-L324):

```
clock_model = estimate_clock_model_with_reroot()
partitions = initialize_partitions()
initialize_marginal(partitions, aln)
run_timetree(no coalescent)                         # first pass
coalescent_tc = optimize_tc()                        # coalescent parameter
run_timetree(with coalescent)                        # second pass

# Iteration loop:
while optimizer.next_iter():
    if coalescent: re-optimize Tc
    run_refinement_iteration():
        relaxed_clock()
        resolve_polytomies()
        update_marginal()                            # ancestral reconstruction
        run_timetree()                               # time inference
    record convergence
```

The refinement iteration at [packages/treetime/src/commands/timetree/refinement.rs#L21-L106](../../packages/treetime/src/commands/timetree/refinement.rs#L21-L106) alternates between `update_marginal()` and `run_timetree()`. Neither `run_optimize_mixed()` nor `initial_guess_mixed()` is called.

This is a known gap: [M-timetree-missing-initial-branch-optimization.md](../port-known-issues/M-timetree-missing-initial-branch-optimization.md). v0 runs `optimize_tree(max_iter=1)` twice before the first time inference. v1 skips this entirely. The infrastructure exists (`PartitionTimetreeOps` extends `PartitionOptimizeOps`) but is not wired.

### Impact of missing ML optimization in timetree

Branch lengths entering the first `run_timetree()` come from the tree file, adjusted by clock regression. For trees from good ML tools (RAxML, IQ-TREE), these branch lengths are already near-optimal and the time inference converges regardless. For trees from parsimony or neighbor-joining, the branch lengths may be poor, causing the time inference to start from a worse state.

v0's `optimize_tree(max_iter=1)` runs one damped round of Brent optimization across all edges before time inference. This is a single-pass refinement, not full convergence. The equivalent in v1 would be:

```rust
initial_guess_mixed(&graph, &mixed_partitions)?;  // optional: set initial from sub counts
run_optimize_mixed(&graph, &mixed_partitions)?;    // one Newton pass over all edges
update_marginal(&graph, &partitions)?;             // recompute profiles with optimized lengths
```

## Advice

### Priority 1: Add outer-loop damping to optimize command

The oscillation known issue ([M-optimize-oscillation-no-damping.md](../port-known-issues/M-optimize-oscillation-no-damping.md)) is the most impactful gap. Two valid approaches:

- **Match v0**: exponential damping `alpha(i) = 1 - 0.75^(i+1)` on branch length updates. Simple, proven for treetime. Requires storing previous branch lengths and blending after each `run_optimize_mixed` call.
- **Match ML tools**: add per-branch likelihood checking (SAFE mode) and step clamping (`max_step = max_bl / max_iters`). More invasive but addresses the root cause (individual edge overshoot) rather than the symptom (oscillation).

v0's approach is simpler and directly addresses the reported problem. Start there.

### Priority 2: Add hard branch length upper bound

v1 is the only tool with no absolute cap on branch lengths. Add `MAX_BRANCH_LENGTH = 4.0` (matching v0) or `10.0` (matching IQ-TREE). Apply as a clamp after Newton/grid search.

### Priority 3: Wire ML optimization into timetree

Add one round of `run_optimize_mixed` before the first `run_timetree()` call, matching v0's `optimize_tree(max_iter=1)`. This seeds time inference with ML-refined branch lengths instead of raw tree-file values. A single undamped pass is sufficient for initial seeding.

### Priority 4: Consider soft Hamming on posteriors (option E) for initial guess

The current hard Hamming (option D) works for well-resolved trees. If testing on datasets with many uncertain internal edges reveals convergence problems (grid search missing optima), replace `edge_subs().len()` with `sum(1 - dot(parent_posterior, child_posterior))` over node posteriors. A third option is the surrogate function `f(t) = a*exp(-b*t) + c` from <a id="cite-3"></a>[Claywell et al., 2017](https://arxiv.org/abs/1706.00659) [[3](#ref-3)], which approximates the likelihood shape for initial bracketing. See [optimize-initial-branch-length-estimation.md](optimize-initial-branch-length-estimation.md) for the full design space.

### Future: Anderson acceleration

<a id="cite-9b"></a>[Henderson & Varadhan, 2018](https://arxiv.org/abs/1803.06673) [[9](#ref-9)] demonstrate that Anderson acceleration with damping and restarts is a drop-in accelerator for EM-like iterations. Applying it to the outer loop could combine the stability of damping with faster convergence. No phylogenetic tool currently uses this approach.

## References

1. <a id="ref-1"></a> Brent RP. "Algorithms for Minimization without Derivatives." Prentice-Hall, 1973. https://maths-people.anu.edu.au/~brent/pd/rpb011i.pdf (accessed 2026-03-24) [↩](#cite-1)
2. <a id="ref-2"></a> Dinh V, Matsen FA IV. "The shape of the one-dimensional phylogenetic likelihood function." 2015. https://arxiv.org/abs/1507.03647 (accessed 2026-03-24) [↩](#cite-2)
3. <a id="ref-3"></a> Claywell BC, Dinh V, McCoy JA, Matsen FA IV. "A surrogate function for one-dimensional phylogenetic likelihoods." 2017. https://arxiv.org/abs/1706.00659 (accessed 2026-03-24) [↩](#cite-3)
4. <a id="ref-4"></a> Fourment M et al. "Automatic differentiation is no panacea for phylogenetic gradient computation." 2022. https://arxiv.org/abs/2211.02168 (accessed 2026-03-24) [↩](#cite-4)
5. <a id="ref-5"></a> Bryant D et al. "LvD: A New Algorithm for Computing the Likelihood of a Phylogeny." 2025. https://arxiv.org/abs/2601.19064 (accessed 2026-03-24) [↩](#cite-5)
6. <a id="ref-6"></a> Clancy D, Lyu H, Roch S. "Sample Complexity of Branch-length Estimation by Maximum Likelihood." 2025. https://arxiv.org/abs/2507.22038 (accessed 2026-03-24) [↩](#cite-6)
7. <a id="ref-7"></a> Mossel E, Roch S, Steel M. "Shrinkage Effect in Ancestral Maximum Likelihood." 2008. https://arxiv.org/abs/0802.0914 (accessed 2026-03-24) [↩](#cite-7)
8. <a id="ref-8"></a> Hobolth A, Yoshida R. "Maximum likelihood estimation of phylogenetic trees via the EM algorithm." 2005. https://arxiv.org/abs/q-bio/0511034 (accessed 2026-03-24) [↩](#cite-8)
9. <a id="ref-9"></a> Henderson NC, Varadhan R. "Damped Anderson Acceleration with Restarts and Applications to Monotone Fixed Point Problems and EM." 2018. https://arxiv.org/abs/1803.06673 (accessed 2026-03-24) [↩¹](#cite-9a) [↩²](#cite-9b)
10. <a id="ref-10"></a> Felsenstein J. "Evolutionary trees from DNA sequences: a maximum likelihood approach." J Mol Evol 17:368-376, 1981. https://doi.org/10.1007/BF01734359 (accessed 2026-03-24) [↩](#cite-10)
11. <a id="ref-11"></a> Pupko T, Pe'er I, Shamir R, Graur D. "A fast algorithm for joint reconstruction of ancestral amino acid sequences." Mol Biol Evol 17:890-896, 2000. https://doi.org/10.1093/oxfordjournals.molbev.a026369 (accessed 2026-03-24) [↩](#cite-11)
12. <a id="ref-12"></a> Stamatakis A. "RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies." Bioinformatics 30(9):1312-1313, 2014. https://doi.org/10.1093/bioinformatics/btu033 (accessed 2026-03-24) [↩](#cite-12)
13. <a id="ref-13"></a> Nguyen LT, Schmidt HA, von Haeseler A, Minh BQ. "IQ-TREE: A fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies." Mol Biol Evol 32(1):268-274, 2015. https://doi.org/10.1093/molbev/msu300 (accessed 2026-03-24) [↩](#cite-13)
14. <a id="ref-14"></a> Minh BQ et al. "IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era." Mol Biol Evol 37(5):1530-1534, 2020. https://doi.org/10.1093/molbev/msaa015 (accessed 2026-03-24) [↩](#cite-14)
15. <a id="ref-15"></a> Sagulenko P, Puller V, Neher RA. "TreeTime: Maximum-likelihood phylodynamic analysis." Virus Evolution 4(1):vex042, 2018. https://doi.org/10.1093/ve/vex042 (accessed 2026-03-24) [↩](#cite-15)
