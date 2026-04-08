# Reliable convergence and method choice for branch length optimization

This document proposes improvements to the `optimize` command's convergence behavior and introduces user-selectable optimization methods. The goal is to match v0's convergence reliability, provide the same systematic method choice already available in the `clock` command, and address three documented known issues.

The proposal is **not accepted** as a whole. P1, P3, P4, and P5 are implemented independently; P2, P6 remain unimplemented. P5 was implemented as part of the `--opt-method` feature (Newton, NewtonSqrt, Brent) with NewtonSqrt as default.

## Problem statement

The `optimize` command has three documented convergence and robustness problems:

1. ~~**Oscillation without damping**~~ (fixed). v1 now applies exponential damping (`--damping`, default 0.75) matching v0.

2. **Initial guess inflated by gaps** (fixed). `initial_guess_mixed()` now filters gap positions from both the substitution count and the effective alignment length denominator.

3. ~~**GTR hardcoded to JC69**~~ (fixed). The `--model` flag is now wired through `get_gtr_sparse()`/`get_gtr_dense()`, matching the ancestral command's dispatch pattern.

Beyond these bugs, the optimize command offers no user control over the per-edge optimization method. The `clock` command exposes Grid, Brent, and Golden Section via `--branch-split-method`. The `optimize` command uses Newton-Raphson with a hardcoded grid search fallback and no way to select alternatives.

## Motivation

### Convergence reliability

The alternating optimization pattern (fix branch lengths and reconstruct ancestral profiles, then fix profiles and optimize branch lengths) is standard in phylogenetic ML inference. Clancy, Lyu & Roch (arXiv 2507.22038) prove exponential convergence in the strong signal-to-noise regime with sufficiently close initialization. Outside that regime, alternation can oscillate: changing all branch lengths shifts marginal profiles, which shifts optimal branch lengths, and so on.

v0 addresses this with exponential damping in `optimize_tree_marginal()` ([packages/legacy/treetime/treetime/treeanc.py#L1297-L1360](../../packages/legacy/treetime/treetime/treeanc.py#L1297-L1360)):

```python
update_val = new_val * (1 - damping ** (i + 1)) + old_val * damping ** (i + 1)
```

This blends new and old branch lengths with a weight that shifts from conservative (75% old at iteration 1) to aggressive (94% new at iteration 10). The damping schedule does not change the fixed point but smooths the path toward it. Polyak & Tremba (arXiv 1703.07810) show that adaptive damping between pure and damped Newton steps gives global convergence guarantees.

v0 also tightens the per-branch Brent tolerance across outer iterations (`tol = 1e-8 + 0.01^(i+1)`), spending less work on early iterations where profiles are still rough.

v1 has none of these safeguards. The convergence loop at [packages/treetime/src/commands/optimize/run.rs#L134-L155](../../packages/treetime/src/commands/optimize/run.rs#L134-L155) applies full Newton-optimal updates each iteration and checks only absolute likelihood change (`--dp`, default 1e-2) for early stopping, which cannot distinguish convergence from oscillation.

### Method choice consistency

The `clock` command provides three optimization methods via `--branch-split-method`:

| Method         | Implementation                                                                                                     |
| -------------- | ------------------------------------------------------------------------------------------------------------------ |
| Grid (default) | `packages/treetime/src/commands/clock/find_best_root/find_best_split.rs`                                           |
| Brent          | `packages/treetime/src/commands/clock/find_best_root/method_brent.rs` using `argmin::BrentOpt`                     |
| Golden Section | `packages/treetime/src/commands/clock/find_best_root/method_golden_section.rs` using `argmin::GoldenSectionSearch` |

The `timetree` command uses Brent for Tc optimization ([packages/treetime/src/commands/timetree/coalescent/optimize_tc.rs#L57-L59](../../packages/treetime/src/commands/timetree/coalescent/optimize_tc.rs#L57-L59)) and for polytomy merge time optimization ([packages/treetime/src/commands/timetree/optimization/polytomy.rs#L243-L245](../../packages/treetime/src/commands/timetree/optimization/polytomy.rs#L243-L245)).

The `optimize` command uses Newton-Raphson with a grid search fallback and no user choice. An end-user cannot request Brent for comparison or debugging, even though the `argmin` dependency and Brent wrapper already exist.

### Practical benefits for end-users

1. **Reliable convergence on diverse datasets.** Damping prevents the optimizer from oscillating on datasets where the marginal profiles and branch lengths are tightly coupled (short alignments, strong compositional bias, many near-zero branches).

2. **Correct results with non-JC69 models.** Wiring `--model` lets users apply appropriate substitution models. JC69 assumes equal base frequencies and equal substitution rates. For GC-biased genomes, amino acid data, or transition/transversion-biased sequences, JC69 produces systematically wrong branch lengths (Hasegawa, Kishino & Yano, 1985).

3. **Debuggability.** When Newton-Raphson produces unexpected results, an end-user can switch to Brent (derivative-free) to check whether the issue is in the derivative computation or the likelihood surface itself.

4. **Consistency across commands.** Users who learn `--branch-split-method` in the clock command expect similar control in the optimize command.

## Proposed changes

### ~~P1: Outer-loop damping~~ (implemented)

Implemented in [packages/treetime/src/commands/optimize/run.rs](../../packages/treetime/src/commands/optimize/run.rs). `--damping` CLI parameter (default 0.75, matching v0). `save_branch_lengths()` captures edge values before `run_optimize_mixed()`, then `apply_damping()` blends `bl_new * (1 - damping^(i+1)) + bl_old * damping^(i+1)`. Value 0.0 disables damping. Value >= 1.0 rejected at parse time.

**Interaction with timetree:** The timetree command has its own iteration loop and does not call `run_optimize()`. If damping proves beneficial, the timetree loop would need independent damping, but that is out of scope for this proposal.

### P2: Progressive tolerance tightening

Add per-iteration tolerance scaling for the Newton convergence criterion. The current Newton inner loop converges when `|delta_bl| <= 0.001 * bl` ([optimize_unified.rs#L258](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L258)). Scale this by iteration:

```
relative_tol = 0.01 * 0.01^i  // tight: 0.01, 0.0001, 0.000001, ...
```

Early iterations use a loose tolerance (profiles are rough, no point in precise branch lengths). Later iterations tighten as profiles stabilize.

**Implementation location:** `packages/treetime/src/commands/optimize/optimize_unified.rs`, pass tolerance as a parameter to `run_optimize_mixed()` (currently at [optimize_unified.rs#L221](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L221)).

### ~~P3: Wire `--model` flag~~ (implemented)

Implemented in [packages/treetime/src/commands/optimize/run.rs](../../packages/treetime/src/commands/optimize/run.rs). Sparse GTR dispatch after `compress_sequences`, dense GTR dispatch after `initialize_marginal` + `update_marginal`, matching the ancestral command's pattern via `get_gtr_sparse()`/`get_gtr_dense()`. A missing `update_marginal` call after dense GTR replacement was added to ensure `initial_guess_mixed` reads fresh edge messages from the real model.

### ~~P4: Fix initial guess gap handling~~ (implemented)

Implemented in `edge_subs()` and `edge_effective_length()` on `PartitionOptimizeOps` ([packages/treetime/src/commands/optimize/partition_ops.rs](../../packages/treetime/src/commands/optimize/partition_ops.rs)). Both sparse and dense implementations filter gap positions. `initial_guess_mixed()` uses per-edge effective alignment length as the denominator.

### P5: Brent as alternative per-edge optimizer

Add Brent's method as an alternative to Newton-Raphson for per-edge branch length optimization. This reuses the existing `argmin::BrentOpt` wrapper from the clock command.

**Why Brent.** Brent's method (Brent, 1973) combines golden section search (guaranteed convergence within bracket) with parabolic interpolation (superlinear convergence, order ~1.325). It needs only function evaluations, not derivatives. This makes it a natural fallback when Newton-Raphson behaves poorly (e.g., near the non-concave region, or when the Hessian is numerically noisy).

**Why not Golden Section separately.** Brent subsumes Golden Section as a component. On the branch length optimization problem - a scalar objective on `[0, MAX_BRANCH_LENGTH]` - Brent dominates Golden Section in both speed and robustness. Golden Section is useful as a standalone method only when the objective lacks the smoothness that parabolic interpolation exploits. The phylogenetic log-likelihood is smooth (infinitely differentiable in `t` for fixed topology), so this case does not arise.

Both RAxML (Stamatakis, 2006; 2014) and IQ-TREE (Nguyen et al., 2015; Minh et al., 2020) use Newton-Raphson as the primary per-branch optimizer with Brent as a fallback, matching the proposed layered approach.

**Proposed CLI interface:**

```
--opt-method newton   (default, current behavior)
--opt-method brent    (derivative-free, guaranteed convergence)
--opt-method grid     (pure grid search, diagnostic)
```

**Implementation approach:**

Create a `BranchLengthOptMethod` enum (Newton, Brent, Grid) analogous to the clock command's `OptimizationMethod`. In `run_optimize_mixed()`:

- **Newton (default):** current code path. Newton-Raphson with clamped step, grid search fallback when Hessian non-negative.
- **Brent:** wrap the likelihood function (without derivatives) in an `argmin` `CostFunction` impl. Bracket: `[0, max(3 * current_bl, 10 * one_mutation)]`. Uses `evaluate_mixed_log_lh_only()` (already exists, computes likelihood without derivatives for efficiency).
- **Grid:** current fallback path, promoted to standalone mode. Evaluate N points on a linear grid and select the maximum. Expose grid size via `--opt-grid-points` (default 100).

**Implementation locations:**

- New args in `packages/treetime/src/commands/optimize/args.rs`
- New enum in `packages/treetime/src/commands/optimize/` (or reuse clock's pattern from `packages/treetime/src/commands/clock/find_best_root/params.rs`)
- Dispatch in `run_optimize_mixed()` at [packages/treetime/src/commands/optimize/optimize_unified.rs#L221-L291](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L221-L291)

### P6: Backtracking line search for Newton

Add an optional backtracking line search within the Newton inner loop. After computing the Newton step `delta = -L'/L''`, verify the step actually improves log-likelihood. If not, halve the step and retry (up to 5 halvings). This prevents individual edge overshoot that can trigger oscillation in the outer loop.

Shea & Schmidt (arXiv 2401.06809) show that Newton with exact line search retains local superlinear convergence while achieving better global progress than backtracking. For per-branch phylogenetic optimization, a simple halving line search is sufficient because the likelihood evaluation is cheap (reuses cached coefficients).

**Implementation location:** [packages/treetime/src/commands/optimize/optimize_unified.rs#L252-L272](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L252-L272), wrap the Newton step in a line search that checks `evaluate_mixed_log_lh_only(contributions, new_bl) > current_log_lh`.

## Proposed feature inventory after implementation

### Per-Edge Optimization

- [x] Newton-Raphson with analytical derivatives (`--opt-method newton`)
- [x] Newton-Raphson in sqrt(t) space (`--opt-method newton-sqrt`, default)
- [x] Brent as alternative method (`--opt-method brent`, P5 implemented)
- [x] Grid search fallback (shared non-concave fallback for Newton methods)
- [ ] Backtracking line search for Newton (P6)
- [x] Zero branch length short-circuit (current)
- [x] Eigenvalue-space coefficient caching (current)

### Convergence Loop

- [x] Iterative reconstruction + optimization (current)
- [x] Early stop on likelihood change (current)
- [x] Exponential damping (P1, implemented)
- [ ] Progressive tolerance tightening (P2)

### GTR Integration

- [x] `--model` wired (P3, implemented)
- [ ] GTR inference in optimization loop (future, not in this proposal)

### Initial Guess

- [x] Gap-aware effective alignment length (P4, implemented)

## Alternative convergence acceleration methods

The outer loop (alternating between marginal reconstruction and branch length optimization) is a fixed-point iteration: `theta_{k+1} = F(theta_k)` where `theta` is the vector of all branch lengths and `F` is "reconstruct profiles, then optimize branch lengths." P1's exponential damping is the simplest acceleration. The optimization literature offers several alternatives with stronger theoretical properties.

### Anderson acceleration (AA)

Anderson acceleration (Anderson, 1965) uses a window of M previous iterates and their residuals to compute the next iterate by solving a constrained least-squares problem. For linear problems it reduces to GMRES (Walker and Ni). Evans et al. (arXiv 1810.08455) prove AA improves the linear convergence rate by a factor related to the residual reduction gain. De Sterck and He (arXiv 2109.14176) show the root-linear convergence factor depends on initial conditions and the acceleration coefficients oscillate near the fixed point.

Henderson and Varadhan (arXiv 1803.06673) propose DAAREM (Damped Anderson Acceleration with Restarts and Monotonicity control), specifically designed for EM-like algorithms. DAAREM combines three safeguards: periodic window restarts to flush stale history, a damping factor blending the AA step with a plain step, and monotonicity enforcement that rejects steps decreasing the objective. DAAREM outperforms SQUAREM on representative EM problems.

Bertrand and Massias (arXiv 2011.10065) apply AA extrapolation directly to coordinate descent, reporting consistent practical speedups over both non-accelerated and Nesterov-accelerated coordinate descent on Lasso and logistic regression.

**Applicability to treetime optimize:** The outer loop is a contractive fixed-point iteration with linear convergence. AA(m) with m=5-10 and monotonicity enforcement (reject steps where total log-likelihood decreases) would be a drop-in replacement for the current damping, with stronger theoretical convergence properties. The per-iteration cost increase is one least-squares solve of size m, negligible compared to the marginal reconstruction. The main implementation cost is maintaining a window of M branch-length vectors.

### SQUAREM

SQUAREM (Varadhan and Roland, 2008; R package: Du and Varadhan, arXiv 1810.11163) accelerates any monotone fixed-point iteration using squared polynomial extrapolation. Each SQUAREM iteration takes two fixed-point steps `theta_1 = F(theta_0)`, `theta_2 = F(theta_1)`, computes `r = theta_1 - theta_0` and `v = (theta_2 - theta_1) - r`, then extrapolates: `theta_sq = theta_0 - 2*alpha*r + alpha^2*v` with step length `alpha = -sqrt(||r||^2 / ||v||^2)`. The objective is checked and if the extrapolated point is worse, the algorithm falls back to `theta_2`.

Typical speedup: 10-50x fewer fixed-point evaluations than plain iteration. SQUAREM uses constant memory (no history window), making it simpler to implement than AA. Three step-length schemes exist (BB-long, BB-short, geometric mean), with the geometric mean (default) performing best.

**Applicability to treetime optimize:** SQUAREM requires two full outer-loop iterations per acceleration step (2x marginal reconstruction + 2x branch length optimization), then one extrapolation + one stabilization evaluation. For the current 20-iteration default, this would mean ~10 SQUAREM steps = ~30 evaluations, vs 20 unaccelerated. The benefit is faster convergence to higher precision, especially for large trees where each outer iteration is expensive.

### Nesterov/momentum acceleration

Classical Nesterov acceleration achieves O(1/k^2) for convex smooth problems. Aujol, Dossal and Rondepierre (arXiv 1805.05719) prove an important caveat: Nesterov acceleration can be _worse_ than plain iteration on sharp (strongly convex) functions because the momentum overshoots the narrow optimum. This is relevant for phylogenetic branch length optimization where likelihood surfaces can have sharp ridges near zero-length branches.

Peng and Yin (arXiv 2405.16020) show that for two-block coordinate descent on least-squares with orthogonal blocks, optimal stepsizes alone (without momentum) achieve convergence rates twice as fast as Polyak momentum. This suggests careful stepsize tuning may be more valuable than momentum for the alternating optimization structure.

**Applicability to treetime optimize:** Nesterov momentum is poorly suited for this problem due to the sharp-function caveat and the non-smooth transitions when branches hit the zero-length boundary. The exponential damping (P1) or Anderson acceleration are safer choices.

### Per-edge Newton safeguards from phylogenetic software

Research on the RAxML-NG and IQ-TREE codebases reveals per-edge Newton safeguards that go beyond v1's current implementation:

**RAxML-NG (via pll-modules):** Six NR variants. The "new NR" (`nr_fast`/`nr_safe`) uses explicit step damping (`dxmax = xmax / max_iters`), bracket tracking with clamping, and concavity safeguards (uses `|H|` when Hessian is non-negative). The `nr_safe` mode recomputes per-branch likelihood after each Newton optimization and reverts the branch length if likelihood worsened. The `FALLBACK` mode starts with `nr_fast` and switches to `nr_safe` when overall likelihood decreases.

**IQ-TREE:** NR implementation numerically identical to PLL "old NR". Uses bisection fallback when Hessian is non-negative or NR step overshoots the bracket. Per-round (not per-branch) monotonicity check: if total likelihood decreases after a full smoothing round, all branch lengths are reverted. If optimized length exceeds 95% of `max_branch_length`, reverts if old length gave better likelihood ("newton raphson diverged, reset").

**v1 comparison:** v1's current NR has step clamping to `[-1.0, current_bl]` and grid search fallback when Hessian is non-negative. It lacks bracket tracking, per-branch likelihood rollback, step damping (`dxmax`), and the bisection fallback. P6's proposed backtracking line search would address the per-branch likelihood check. Bracket tracking and step damping could be added as independent improvements.

### Comparison summary

| Method                            | Memory              | Per-step cost                     | Monotonicity | Convergence rate improvement    | Complexity |
| --------------------------------- | ------------------- | --------------------------------- | ------------ | ------------------------------- | ---------- |
| Exponential damping (P1, current) | O(n) branch lengths | 1 outer iteration                 | No guarantee | Smooths path, same fixed point  | Minimal    |
| Backtracking line search (P6)     | O(1) per edge       | 1-5 likelihood evals per edge     | Per-edge     | Prevents overshoot              | Low        |
| Progressive tolerance (P2)        | O(1)                | Fewer inner NR iters early        | No           | Reduces wasted work             | Minimal    |
| SQUAREM                           | O(n) branch lengths | 2-3 outer iterations              | Configurable | 10-50x fewer total iterations   | Low        |
| Anderson acceleration (DAAREM)    | O(mn) for window m  | 1 outer iteration + O(nm^2) solve | Enforced     | Provably improved linear rate   | Medium     |
| Nesterov momentum                 | O(n) branch lengths | 1 outer iteration                 | No guarantee | O(1/k^2) convex, worse on sharp | Low        |

n = number of branches, m = AA window size (5-10 typical).

## Priority ordering

| Priority | Item | Rationale                                                    |
| -------- | ---- | ------------------------------------------------------------ |
| 1        | P1   | Directly fixes the reported oscillation problem              |
| 2        | P3   | Consistency with ancestral command, infrastructure exists    |
| 3        | P4   | Fixes incorrect initial guess for gappy alignments           |
| 4        | P5   | Method choice consistency with clock command                 |
| 5        | P2   | Refinement of P1, reduces wasted work on early iterations    |
| 6        | P6   | Refinement of Newton step, diminishing returns if P1 is done |

P1 and P3 can be implemented independently. P4 is independent of all others. P5 requires no changes to the Newton path (additive). P2 and P6 refine the Newton path and interact with P1 (damping + tolerance + line search together control step behavior).

SQUAREM or DAAREM could replace P1 in the future for faster convergence on large trees. Both are drop-in replacements that operate on the same outer loop structure. SQUAREM is simpler to implement. DAAREM has stronger theoretical backing and monotonicity guarantees.

## Validation plan

1. **Oscillation test.** Run optimize on datasets known to oscillate (ebola, flu/h3n2/200). Compare likelihood traces with and without damping. Damped runs should show monotonic or near-monotonic likelihood improvement.

2. **Model parity test.** Run optimize with `--model=hky85` on datasets with known transition/transversion bias. Compare branch lengths against ancestral command output with the same model.

3. **Gap initialization test.** Compare initial branch length estimates on gappy alignments (tb, lassa/L) before and after P4. Verify estimates are closer to v0's gap-excluded Hamming distance.

4. **Method equivalence test.** Run optimize with `--opt-method=newton` and `--opt-method=brent` on the same inputs. Final branch lengths should agree within tolerance (1e-4 relative difference). Likelihood should be equal within 1e-6.

5. **Regression test.** Run smoke tests on all datasets with default parameters. Verify no regression in likelihood or branch length output compared to pre-change baseline.

## References

- Anderson DG (1965). Iterative procedures for nonlinear integral equations. J ACM 12(4):547-560.
- Aujol JF, Dossal C, Rondepierre A (2019). Optimal convergence rates for Nesterov acceleration. arXiv 1805.05719.
- Bertrand Q, Massias M (2020). Anderson acceleration of coordinate descent. arXiv 2011.10065.
- Brent RP (1973). _Algorithms for Minimization without Derivatives_. Prentice-Hall.
- Clancy D, Lyu S, Roch S (2025). Sample complexity of branch-length estimation by maximum likelihood. arXiv 2507.22038.
- De Sterck H, He Y (2021). Linear asymptotic convergence of Anderson acceleration: fixed-point analysis. arXiv 2109.14176.
- Dinh V, Matsen FA IV (2015). The shape of the one-dimensional phylogenetic likelihood function. arXiv 1507.03647.
- Du Y, Varadhan R (2018). SQUAREM: An R package for off-the-shelf acceleration of EM, MM and other EM-like monotone algorithms. arXiv 1810.11163.
- Evans C, Pollock S, Rebholz LG, Xiao M (2018). A proof that Anderson acceleration improves the convergence rate in linearly converging fixed-point methods. arXiv 1810.08455.
- Hasegawa M, Kishino H, Yano T (1985). Dating of the human-ape splitting by a molecular clock of mitochondrial DNA. J Mol Evol 22:160-174.
- Henderson NC, Varadhan R (2019). Damped Anderson acceleration with restarts and monotonicity control for accelerating EM and EM-like algorithms. arXiv 1803.06673.
- Kenney T, Gu H (2012). Hessian calculation for phylogenetic likelihood based on the pruning algorithm. Statistical Applications in Genetics and Molecular Biology 11(4):14.
- Kozlov AM, Darriba D, Flouri T, Morel B, Stamatakis A (2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics 35(21):4453-4455.
- Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Mol Biol Evol 37(5):1530-1534.
- Nguyen LT, Schmidt HA, von Haeseler A, Minh BQ (2015). IQ-TREE: A fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Mol Biol Evol 32(1):268-274.
- Peng W, Yin W (2024). Optimal stepsizes for two-block alternating gradient descent. arXiv 2405.16020.
- Polyak BT, Tremba AA (2017). New versions of Newton method: step-size choice, convergence domain and under-determined equations. arXiv 1703.07810.
- Shea D, Schmidt M (2024). Greedy Newton: Newton's method with exact line search. arXiv 2401.06809.
- Stamatakis A (2006). RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models. Bioinformatics 22(21):2688-2690.
- Stamatakis A (2014). RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics 30(9):1312-1313.
- Varadhan R, Roland C (2008). Simple and globally convergent methods for accelerating the convergence of any EM algorithm. Scandinavian Journal of Statistics 35(2):335-353.
