# Branch-length optimization: convergence and robustness

## Problem statement

The `optimize` command alternates marginal ancestral reconstruction (E-step) and per-edge branch-length optimization (M-step) in a loop. This is a form of ECM (Expectation Conditional Maximization; Meng and Rubin 1993), standard in phylogenetic ML inference (Felsenstein 1981).

Dense mode operates as true soft-EM: every position stores a full probability vector, the E-step is purely continuous, and the objective is guaranteed non-decreasing per iteration (Dempster, Laird, and Rubin 1977; Wu 1983). Dense mode converges reliably.

Sparse mode introduces an approximation: positions are classified as "variable" (individual posterior) or "fixed" (canonical-state profile bucketed by multiplicity). This classification can change between iterations, producing a discrete jump in the per-edge objective. The sparse loop is therefore not a standard EM and does not have the monotone convergence guarantee. The immediate non-convergence bug is addressed in M-optimize-sparse-em-2-cycle (resolved) by making the convergence loop robust to non-monotone behavior (three-condition convergence check, damping floor, v0 defaults).

This proposal addresses deeper architectural improvements that would eliminate or reduce the sparse approximation error, improve convergence speed, and add per-edge robustness features present in other phylogenetic software.

## Architecture overview

### Outer loop

`run_optimize_loop()` [packages/treetime/src/commands/optimize/run.rs#L275](../../packages/treetime/src/commands/optimize/run.rs#L275)

```
for i in 0..max_iter:
  1. E-step: update_marginal(sparse) + update_marginal(dense)
  2. convergence check
  3. M-step: run_optimize_mixed (Newton/Brent per edge, with indel Poisson term)
  4. find zero-optimal internal edges (before damping)
  5. apply_damping
  6. topology cleanup: collapse zero edges, merge shared mutations
```

### Sparse representation

The sparse E-step ([packages/treetime/src/representation/partition/marginal_passes.rs](../../packages/treetime/src/representation/partition/marginal_passes.rs)) propagates messages carrying:

- `variable: BTreeMap<usize, VarPos>` -- per-position posterior `dis: Array1<f64>` and reference state
- `fixed: BTreeMap<AsciiChar, Array1<f64>>` -- one propagated profile per canonical state, shared by all same-state positions
- `fixed_counts: Composition` -- multiplicity per canonical state

The variable set is rebuilt from scratch each iteration and can grow or shrink.

The sparse M-step evaluator ([packages/treetime/src/commands/optimize/optimize_sparse.rs](../../packages/treetime/src/commands/optimize/optimize_sparse.rs)) computes eigendecomposition coefficients `parent.dot(V) * child.dot(V_inv.T)` from individual posteriors for variable positions and from bucketed canonical profiles for fixed positions. The coefficient computation uses soft distributions (not hard MAP states), but the variable/fixed classification boundary produces a discontinuity: the same position yields different coefficients depending on which path evaluates it.

### Dense representation

The dense E-step ([packages/treetime/src/representation/partition/marginal_dense.rs](../../packages/treetime/src/representation/partition/marginal_dense.rs)) stores full `Array2<f64>` probability matrices. The dense evaluator ([packages/treetime/src/commands/optimize/optimize_dense.rs](../../packages/treetime/src/commands/optimize/optimize_dense.rs)) computes coefficients from these matrices with multiplicity 1 per position. No variable/fixed classification. True soft-EM.

### Indel contribution

[packages/treetime/src/commands/optimize/optimize_indel.rs](../../packages/treetime/src/commands/optimize/optimize_indel.rs)

Per-edge Poisson term $\ell_{\text{indel}}(t) = k \ln(\mu t) - \mu t - \ln(k!)$ where $t$ is branch length, $k$ is the observed indel count on the edge, and $\mu$ is the global indel rate. The global rate $\hat{\mu} = \sum_e k_e / \sum_e t_e$ is recomputed each iteration, where the sums run over all edges $e$. v1-only feature (v0 ignores indels). See [optimize-indel-contribution-to-likelihood](../port-intentional-changes/optimize-indel-contribution-to-likelihood.md).

### Per-edge optimizers

Six methods via `--opt-method`: Newton and Brent in $t$, $\sqrt{t}$, $\ln(t)$ spaces. Default `brent-sqrt` (matches v0). Grid search fallback for non-concave regions. See [optimize-newton-raphson-per-edge](../port-intentional-changes/optimize-newton-raphson-per-edge.md).

## Implemented features

| ID  | Feature                                               | Reference                                                                                                                                                                                                        |
| --- | ----------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| I1  | Exponential damping (`--damping`, default 0.75)       | `apply_damping()` [commands/optimize/run.rs#L513](../../packages/treetime/src/commands/optimize/run.rs#L513)                                                                                                     |
| I2  | Six per-edge methods (`--opt-method`)                 | [commands/optimize/method_newton.rs](../../packages/treetime/src/commands/optimize/method_newton.rs), [commands/optimize/method_brent.rs](../../packages/treetime/src/commands/optimize/method_brent.rs)         |
| I3  | Grid search fallback (100-point, non-concave regions) | [commands/optimize/optimize_unified.rs](../../packages/treetime/src/commands/optimize/optimize_unified.rs)                                                                                                       |
| I4  | `--model` wired to GTR dispatch                       | [commands/optimize/run.rs](../../packages/treetime/src/commands/optimize/run.rs)                                                                                                                                 |
| I5  | Gap-aware initial guess (`edge_effective_length`)     | [commands/optimize/partition_ops.rs](../../packages/treetime/src/commands/optimize/partition_ops.rs)                                                                                                             |
| I6  | Zero-branch short-circuit (unimodal models)           | `is_zero_branch_optimal()` [commands/optimize/optimize_unified.rs#L328](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L328)                                                                  |
| I7  | Eigenvalue-space coefficient caching                  | [commands/optimize/optimize_sparse.rs](../../packages/treetime/src/commands/optimize/optimize_sparse.rs), [commands/optimize/optimize_dense.rs](../../packages/treetime/src/commands/optimize/optimize_dense.rs) |
| I8  | Poisson indel contribution                            | [commands/optimize/optimize_indel.rs](../../packages/treetime/src/commands/optimize/optimize_indel.rs)                                                                                                           |
| I9  | Collapse guard (skip edges with mutations)            | `find_zero_optimal_internal_edges()` [commands/optimize/run.rs#L563](../../packages/treetime/src/commands/optimize/run.rs#L563)                                                                                  |

## Proposed improvements

### Sparse convergence

**P1. Monotonic variable position set.** Once a position enters the variable set on any edge, it remains variable in subsequent iterations. The distributions at those positions continue updating; only set membership is frozen.

- Addresses: the discrete variable/fixed reclassification that causes the likelihood 2-cycle
- Evidence: variable count oscillates 35095/35109 on sc2/2844 in lockstep with the 2-cycle; stable at 1121 on flu/h3n2/20 where convergence is normal
- Implementation: ~30 lines in [representation/partition/marginal_passes.rs](../../packages/treetime/src/representation/partition/marginal_passes.rs). Store previous variable keys per edge. Reinstate previously-variable positions with the canonical-state profile if absent from the new variable map.
- Pros: eliminates the discrete jump at its source. Minimal performance cost (~14 extra positions out of ~35000).
- Cons: makes the sparse approximation history-dependent -- which positions are treated individually depends on early iterations, not only on the current posterior. The ~14 oscillating positions are genuinely ambiguous (near-equal posterior for two states), so promoting them to individual evaluation is the more accurate classification. After 1-2 iterations the set stabilizes.
- Note: this is a model change, not a neutral fix. The immediate bug fix in M-optimize-sparse-em-2-cycle (resolved) addresses convergence through loop robustness instead, without changing the sparse model.

**P2. Transparent variable/fixed boundary.** Ensure that a position's eigendecomposition coefficient is numerically identical regardless of whether it is evaluated via the variable path (individual posterior) or the fixed path (bucketed canonical profile).

- Addresses: the fundamental source of the discontinuity. If both paths produce the same coefficient, reclassification is harmless.
- Approach: for fixed positions, verify that the bucketed canonical profile matches what the individual posterior would produce. Positions where they diverge are promoted to variable.
- Pros: the variable/fixed boundary becomes a pure performance optimization with no numerical effect. No history-dependence concern.
- Cons: requires determining which fixed positions are "truly invariant" (canonical profile equals individual posterior) versus "approximately invariant" (profiles differ due to nearby variation). Nontrivial analysis of the sparse message structure.
- Depends on: detailed analysis of when the two paths diverge numerically. The fixed profile is a propagated canonical state vector shared by all same-state positions; the variable profile is position-specific posterior incorporating the full tree context. These are inherently different quantities for positions adjacent to variable sites.

**P3. Per-position sparse evaluator.** Eliminate the fixed bucket entirely in the M-step. Evaluate every position individually using its actual posterior distribution.

- Addresses: the variable/fixed boundary by removing it
- Approach: variable positions use individual posteriors (current). Invariant positions use canonical-state profiles computed from the alphabet (fast, one dot product per position). No bucketing by state, no multiplicity.
- Pros: makes sparse M-step evaluation equivalent to dense for correctness. Eliminates the boundary entirely.
- Cons: evaluates all ~30,000 positions instead of ~200 variable + ~5 fixed buckets. Each position is one dot product (fast), but ~150x more positions than current sparse. Sparse-over-dense speedup would shrink from ~150x to ~2-5x.
- Alternative: apply per-position evaluation only in the M-step; keep bucketed propagation in the E-step where the speedup is largest.

### Outer loop

**P4. Indel rate caching.** Compute `estimate_indel_rate` once (after topology stabilizes) and pass to `run_optimize_mixed` as a fixed parameter.

- Addresses: the $\hat{\mu}$ feedback loop that amplifies the 2-cycle. On sc2/2844, $\hat{\mu} \approx 12{,}000$ oscillates proportionally with total BL.
- Implementation: ~10 lines. Add `indel_rate: f64` parameter to `run_optimize_mixed`. Recompute in `run_optimize_loop` while topology changes occur (`prune_and_merge_in_loop` returns true); cache once topology stabilizes.
- Pros: eliminates feedback between total BL and indel rate in the stable regime
- Cons: changes the optimization objective. Under the current model, $\hat{\mu}$ is defined from current branch lengths and varies between iterations. Caching freezes it, optimizing under a different (fixed-rate) model. In practice the rate converges along with branch lengths, so the fixed-rate model is accurate at convergence. But the path to convergence differs.
- Note: this is a model change, not neutral hardening.
- Rejected storage location: storing the indel rate on the GTR struct (suggested as a way to make it visible during marginal passes) conflates the substitution model with the indel counting process. The GTR struct parameterizes a Markov chain (rate matrix, equilibrium frequencies, eigendecomposition); the indel rate parameterizes an independent Poisson process on branch lengths. A partition-level field or dedicated `IndelModel` struct is preferable if persistence is needed.

**P5. SQUAREM acceleration.** Replace exponential damping with squared polynomial extrapolation (Varadhan and Roland 2008).

Each SQUAREM step takes two plain EM steps $\theta_1 = F(\theta_0)$, $\theta_2 = F(\theta_1)$, computes $r = \theta_1 - \theta_0$ and $v = (\theta_2 - \theta_1) - r$, then extrapolates:

$$\theta_{\text{sq}} = \theta_0 - 2\alpha\, r + \alpha^2 v$$

where $\theta_i$ are parameter vectors at successive EM iterations, $F$ is the EM fixed-point map, $r$ and $v$ are first- and second-order displacement vectors, and $\alpha = -\sqrt{\|r\|^2 / \|v\|^2}$ is the step-length ratio. Falls back to $\theta_2$ if the extrapolated point worsens the objective.

- Typical speedup: 10-50x fewer fixed-point evaluations. Constant memory.
- Pros: drop-in replacement for damping. Simple implementation (~50 lines). Provably convergent for monotone maps. Built-in fallback on non-monotonicity.
- Cons: two E+M evaluations per acceleration step. Requires P1 or P3 for the sparse loop to be approximately monotone.
- Depends on: P1 or P3 (SQUAREM assumes monotone maps; the fallback handles small violations but large discrete jumps can destabilize the step-length estimate)

**P6. DAAREM acceleration.** Anderson acceleration with restarts and monotonicity enforcement (Henderson and Varadhan 2019).

Uses a window of $m$ previous iterates ($m = 5\text{--}10$) and residuals to solve a constrained least-squares problem, where $m$ is the window size. Combines periodic restarts, damping, and monotonicity enforcement.

- Pros: provably improved linear convergence rate (Evans et al. 2020). Monotonicity enforcement handles non-monotone objectives directly, making it more robust than SQUAREM for the sparse loop.
- Cons: window memory $O(mn)$ where $n$ is the parameter dimension, more complex implementation (~100 lines + least-squares solver).
- Priority: lower than P5 (SQUAREM is simpler and sufficient for typical datasets)

### Per-edge optimization

**P7. Progressive tolerance tightening.** Scale per-branch optimizer tolerance by outer iteration: $\text{tol} = \text{tol}_{\text{base}} \cdot 0.01^{i}$ where $i$ is the iteration index and $\text{tol}_{\text{base}}$ is the base tolerance.

v0 uses $\text{tol} = 10^{-8} + 0.01^{i+1}$ in `optimal_t_compressed` ([packages/legacy/treetime/treetime/gtr.py#L816-L920](../../packages/legacy/treetime/treetime/gtr.py#L816-L920)). Early iterations use loose tolerance (profiles are rough). Later iterations tighten as profiles stabilize.

- Implementation: add `tolerance` parameter to `run_optimize_mixed`, pass from `run_optimize_loop` with iteration-dependent scaling
- Pros: reduces wasted computation on early iterations. Matches v0 behavior.
- Cons: minimal impact for typical datasets (loop converges in 5-10 iterations)

**P8. Backtracking line search for Newton.** After computing the Newton step $\Delta = -f'/f''$ where $f'$ and $f''$ are the first and second derivatives of the log-likelihood with respect to branch length, verify the step improves log-likelihood. Halve and retry (up to 5 halvings) if not.

RAxML-NG's `nr_safe` mode (Kozlov et al. 2019) recomputes per-branch likelihood and reverts if worsened. IQ-TREE (Minh et al. 2020) checks per-round monotonicity. v1 has neither.

- Implementation: wrap the Newton step with a likelihood-only evaluation check
- Pros: prevents per-edge overshoot. Matches RAxML-NG safeguard.
- Cons: 1-5 extra likelihood evaluations per edge per iteration (small cost, reuses cached coefficients)

### Model integration

**P9. GTR inference in optimization loop.** Re-estimate the GTR model per iteration, matching v0's `infer_gtr` parameter.

v0 code: [packages/legacy/treetime/treetime/treeanc.py#L1307-L1309](../../packages/legacy/treetime/treetime/treeanc.py#L1307-L1309).

- Implementation: add `--infer-gtr` flag, call `infer_gtr_sparse`/`infer_gtr_dense` + `update_marginal` per iteration
- Pros: joint optimization of model parameters and branch lengths. Required for `--model infer` in the optimize command.
- Cons: adds a third block to the ECM iteration. Preserves convergence guarantees (Meng and Rubin 1993) but may slow convergence rate.
- Related: [M-gtr-sparse-composition-stale-after-marginal](../port-known-issues/M-gtr-sparse-composition-stale-after-marginal.md)

## Rejected approaches

| Approach                            | Reason                                                                                                                                                                                                                                                                                                                                                                                                            |
| ----------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Nesterov/momentum acceleration      | Aujol, Dossal, and Rondepierre (2019): can be worse than plain iteration on sharp functions. Phylogenetic BL surfaces have sharp ridges near $t = 0$. Peng and Yin (2024): optimal stepsizes outperform momentum for two-block coordinate descent.                                                                                                                                                                |
| Golden Section as standalone method | Brent subsumes it. No scenario where Golden Section alone is preferable on the smooth phylogenetic likelihood.                                                                                                                                                                                                                                                                                                    |
| Copy v0's signed convergence check  | v0 uses `deltaLH < LHtol` (signed), accepting any likelihood decrease as convergence. This conflates convergence with failure and works in v0 only because dense soft-EM with damping produces near-monotone progress. See [optimize-signed-convergence-check](../port-v0-errata/optimize-signed-convergence-check.md). v1's three-condition check (converged, oscillating, worsened-with-revert) is more robust. |

## Priority

| Tier                | Items      | Rationale                                                                                                   |
| ------------------- | ---------- | ----------------------------------------------------------------------------------------------------------- |
| Sparse correctness  | P1, P2, P3 | Reduce or eliminate the variable/fixed boundary error. P1 is simplest, P3 eliminates the boundary entirely. |
| Convergence speed   | P5, P6     | Faster convergence for large trees. SQUAREM is simpler.                                                     |
| v0 parity           | P7, P9     | Match v0 features (progressive tolerance, GTR inference in loop).                                           |
| Per-edge robustness | P8         | Match RAxML-NG/IQ-TREE safeguards.                                                                          |
| Indel model         | P4         | Remove feedback amplifier. Model change, requires careful validation.                                       |

## Comparison with other phylogenetic software

| Feature              | v0           | v1 (current)                               | RAxML-NG                     | IQ-TREE            |
| -------------------- | ------------ | ------------------------------------------ | ---------------------------- | ------------------ |
| Per-edge method      | Brent (sqrt) | 6 methods                                  | Newton (6 variants)          | Newton             |
| Outer damping        | exponential  | exponential + floor                        | none (per-edge step damping) | none               |
| Convergence check    | signed delta | abs delta + oscillation detection + revert | per-round rollback           | per-round rollback |
| Tolerance tightening | yes          | no                                         | no                           | no                 |
| Per-branch rollback  | no           | no                                         | yes (`nr_safe`)              | no                 |
| Indel contribution   | no           | yes (Poisson)                              | no                           | no                 |
| GTR re-estimation    | yes          | no                                         | yes                          | yes                |
| Sparse/dense modes   | dense only   | both                                       | dense only                   | dense only         |

## Validation

- **Dense/sparse equivalence**: compare optimized branch lengths between `--dense true` and `--dense false` on flu/h3n2/20 and sc2/2844. After P2 or P3, results should agree within floating-point tolerance.
- **Convergence on sc2/2844**: after P1 or P3, the variable position count should be non-decreasing across iterations and the likelihood should be approximately monotone.
- **v0 parity**: compare v1 output with v0 on sc2/2844 (golden master). See [M-optimize-gm-per-branch-divergence](../port-known-issues/M-optimize-gm-per-branch-divergence.md).
- **No regression**: flu/h3n2/20 should converge identically before and after changes.
- **Acceleration**: with P5, convergence to the same optimum in fewer iterations on ebola, flu/h3n2/200.

## Historical context

The optimizer convergence work progressed through three phases:

1. **Initial convergence fixes**: damping (I1), gap-aware initial guess (I5), `--model` wiring (I4), six per-edge methods (I2). Damping alone is insufficient (see phase 2).

2. **Sparse 2-cycle discovery**: investigation of non-convergence on sc2/2844 (M-optimize-sparse-em-2-cycle (resolved)) revealed that damping alone is insufficient. The sparse variable/fixed position classification oscillates between iterations, creating a discrete objective discontinuity. Dense mode has no such boundary and converges correctly. PR [neherlab/treetime#558](https://github.com/neherlab/treetime/pull/558) commented out `estimate_indel_rate` in `initial_guess_mixed` based on incorrect root-cause analysis. v1 defaults at the time (`max_iter=20`, `dp=0.01`) diverged from v0 (`max_iter=10`, `LHtol=0.1`) without documented reason (now aligned). Peer review identified that v0's signed convergence check is a defect (see [errata](../port-v0-errata/optimize-signed-convergence-check.md)) and that variable-set freezing and indel-rate caching are model changes requiring separate validation.

3. **Immediate fix** (M-optimize-sparse-em-2-cycle (resolved)): makes the convergence loop robust to non-monotone behavior without changing the sparse model. This proposal covers the architectural improvements that address the non-monotonicity at its source.

## Related

- M-optimize-sparse-em-2-cycle (resolved) -- the immediate bug fix
- [M-optimize-gm-per-branch-divergence](../port-known-issues/M-optimize-gm-per-branch-divergence.md) -- per-branch v0 parity
- [M-gtr-sparse-composition-stale-after-marginal](../port-known-issues/M-gtr-sparse-composition-stale-after-marginal.md) -- stale Fitch composition
- [optimize-signed-convergence-check](../port-v0-errata/optimize-signed-convergence-check.md) -- v0 erratum: signed convergence check
- [optimize-indel-model-alternatives](optimize-indel-model-alternatives.md) -- alternative indel models (orthogonal)
- [optimize-indel-contribution-to-likelihood](../port-intentional-changes/optimize-indel-contribution-to-likelihood.md) -- v1-only Poisson indel term
- [optimize-newton-raphson-per-edge](../port-intentional-changes/optimize-newton-raphson-per-edge.md) -- per-method rationale
- [docs/reports/iterative-tree-refinement/](../reports/iterative-tree-refinement/_index.md) -- damping, convergence theory, loop variants
- [docs/algorithms/optimize.md](../algorithms/optimize.md) -- design document

## References

- Aujol, J.-F., C. Dossal, and A. Rondepierre. 2019. "Optimal Convergence Rates for Nesterov Acceleration." SIAM Journal on Optimization 29(4):3131-3153. https://doi.org/10.1137/18m1186757
- Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. "Maximum Likelihood from Incomplete Data via the EM Algorithm." Journal of the Royal Statistical Society Series B: Statistical Methodology 39(1):1-22. https://doi.org/10.1111/j.2517-6161.1977.tb01600.x
- Evans, C., S. Pollock, L. G. Rebholz, and M. Xiao. 2020. "A Proof That Anderson Acceleration Improves the Convergence Rate in Linearly Converging Fixed-Point Methods (But Not in Those Converging Quadratically)." SIAM Journal on Numerical Analysis 58(1):788-810. https://doi.org/10.1137/19m1245384
- Felsenstein, J. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." Journal of Molecular Evolution 17(6):368-376. https://doi.org/10.1007/BF01734359
- Henderson, N. C., and R. Varadhan. 2019. "Damped Anderson Acceleration With Restarts and Monotonicity Control for Accelerating EM and EM-like Algorithms." Journal of Computational and Graphical Statistics 28(4):834-846. https://doi.org/10.1080/10618600.2019.1594835
- Kozlov, A. M., et al. 2019. "RAxML-NG: A Fast, Scalable and User-Friendly Tool for Maximum Likelihood Phylogenetic Inference." Bioinformatics 35(21):4453-4455. https://doi.org/10.1093/bioinformatics/btz305
- Meng, X.-L., and D. B. Rubin. 1993. "Maximum Likelihood Estimation via the ECM Algorithm: A General Framework." Biometrika 80(2):267-278. https://doi.org/10.1093/biomet/80.2.267
- Minh, B. Q., et al. 2020. "IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era." Molecular Biology and Evolution 37(5):1530-1534. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7182206/ (DOI: https://doi.org/10.1093/molbev/msaa015)
- Peng, W., and W. Yin. 2024. "Optimal Stepsizes for Two-Block Alternating Gradient Descent." arXiv 2405.16020. https://arxiv.org/abs/2405.16020
- Varadhan, R., and C. Roland. 2008. "Simple and Globally Convergent Methods for Accelerating the Convergence of Any EM Algorithm." Scandinavian Journal of Statistics 35(2):335-353. https://doi.org/10.1111/j.1467-9469.2007.00585.x
- Wu, C. F. Jeff. 1983. "On the Convergence Properties of the EM Algorithm." Annals of Statistics 11(1):95-103. https://projecteuclid.org/journals/annals-of-statistics/volume-11/issue-1/On-the-Convergence-Properties-of-the-EM-Algorithm/10.1214/aos/1176346060.full (DOI: https://doi.org/10.1214/aos/1176346060)
