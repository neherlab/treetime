# Reliable convergence and method choice for branch length optimization

This document proposes improvements to the `optimize` command's convergence behavior and introduces user-selectable optimization methods. The goal is to match v0's convergence reliability, provide the same systematic method choice already available in the `clock` command, and address three documented known issues.

The proposal is **not accepted** and **not implemented**.

## Problem statement

The `optimize` command has three documented convergence and robustness problems:

1. **Oscillation without damping** ([M-optimize-oscillation-no-damping](../port-known-issues/M-optimize-oscillation-no-damping.md)). The outer loop alternates between marginal reconstruction and per-edge optimization without blending new and old branch lengths. v0 applies exponential damping (`damping=0.75`); v1 does not.

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

v1 has none of these safeguards. The convergence loop at [packages/treetime/src/commands/optimize/run.rs#L121-L141](../../packages/treetime/src/commands/optimize/run.rs#L121-L141) applies full Newton-optimal updates each iteration and checks only absolute likelihood change (`--dp`, default 1e-2) for early stopping, which cannot distinguish convergence from oscillation.

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

### P1: Outer-loop damping

Add exponential damping to the convergence loop. After `run_optimize_mixed()` updates each edge's branch length, blend with the previous value:

```
bl_new = bl_optimized * (1 - damping^(i+1)) + bl_old * damping^(i+1)
```

Expose `--damping` CLI parameter (default 0.75, matching v0). Value 0.0 disables damping (current behavior). Value 1.0 would prevent any update (degenerate, reject at parse time).

**Implementation location:** `packages/treetime/src/commands/optimize/run.rs`, inside the convergence loop at lines 121-141. Save branch lengths before `run_optimize_mixed()`, then blend after.

**Interaction with timetree:** The timetree command has its own iteration loop and does not call `run_optimize()`. If damping proves beneficial, the timetree loop would need independent damping, but that is out of scope for this proposal.

### P2: Progressive tolerance tightening

Add per-iteration tolerance scaling for the Newton convergence criterion. The current Newton inner loop converges when `|delta_bl| <= 0.001 * bl` ([optimize_unified.rs#L217](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L217)). Scale this by iteration:

```
relative_tol = 0.01 * 0.01^i  // tight: 0.01, 0.0001, 0.000001, ...
```

Early iterations use a loose tolerance (profiles are rough, no point in precise branch lengths). Later iterations tighten as profiles stabilize.

**Implementation location:** `packages/treetime/src/commands/optimize/optimize_unified.rs`, pass tolerance as a parameter to `run_optimize_mixed()`.

### P3: Wire `--model` flag

Follow the ancestral command's pattern ([packages/treetime/src/commands/ancestral/run.rs#L128-L178](../../packages/treetime/src/commands/ancestral/run.rs#L128-L178)):

1. Create partitions with dummy JC69 (existing code, no change)
2. Run Fitch compression / `initialize_marginal()` (existing code, no change)
3. Replace dummy GTR with `get_gtr_sparse(model_name, partition, &graph)` or `get_gtr_dense(model_name, partition, &graph)` (new: wire `model_name` arg)
4. For `--model=infer` on dense: run preliminary marginal pass with dummy GTR to populate profiles, infer real GTR, then re-run (matching ancestral command's two-pass approach)

All infrastructure exists: `get_gtr_by_name()`, `get_gtr_sparse()`, `get_gtr_dense()`, `infer_gtr_sparse()`, `infer_gtr_dense()` in `packages/treetime/src/gtr/get_gtr.rs`.

**Implementation location:** `packages/treetime/src/commands/optimize/run.rs`, replace the four `jc69(JC69Params::default())?` calls at lines 65, 81, 94, 105 with model dispatch calls.

### P4: Fix initial guess gap handling

In `initial_guess_mixed()` ([packages/treetime/src/commands/optimize/optimize_unified.rs#L257-L280](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L257-L280)), compute per-edge effective alignment length excluding positions where either parent or child has a gap or ambiguous state. Use this as the denominator instead of the full partition length.

For sparse partitions, filter `edge_subs()` to exclude substitutions where either the reference or query state maps to a gap character. For dense partitions, compute effective length from the per-edge message intersection.

v0's `state_pair()` ([packages/legacy/treetime/treetime/gtr.py#L631-L712](../../packages/legacy/treetime/treetime/gtr.py#L631-L712)) accepts `ignore_gaps=True` (the default), which excludes gap positions from both the substitution count and the total pair count.

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
- Dispatch in `run_optimize_mixed()` at [packages/treetime/src/commands/optimize/optimize_unified.rs#L180-L250](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L180-L250)

### P6: Backtracking line search for Newton

Add an optional backtracking line search within the Newton inner loop. After computing the Newton step `delta = -L'/L''`, verify the step actually improves log-likelihood. If not, halve the step and retry (up to 5 halvings). This prevents individual edge overshoot that can trigger oscillation in the outer loop.

Shea & Schmidt (arXiv 2401.06809) show that Newton with exact line search retains local superlinear convergence while achieving better global progress than backtracking. For per-branch phylogenetic optimization, a simple halving line search is sufficient because the likelihood evaluation is cheap (reuses cached coefficients).

**Implementation location:** [packages/treetime/src/commands/optimize/optimize_unified.rs#L211-L231](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L211-L231), wrap the Newton step in a line search that checks `evaluate_mixed_log_lh_only(contributions, new_bl) > current_log_lh`.

## Proposed feature inventory after implementation

### Per-Edge Optimization

- [x] Newton-Raphson with analytical derivatives (current, P5 keeps as default)
- [x] Grid search fallback (current, P5 promotes to standalone mode)
- [ ] Brent as alternative method (P5)
- [ ] Backtracking line search for Newton (P6)
- [x] Zero branch length short-circuit (current)
- [x] Eigenvalue-space coefficient caching (current)

### Convergence Loop

- [x] Iterative reconstruction + optimization (current)
- [x] Early stop on likelihood change (current)
- [ ] Exponential damping (P1)
- [ ] Progressive tolerance tightening (P2)

### GTR Integration

- [ ] `--model` wired (P3)
- [ ] GTR inference in optimization loop (future, not in this proposal)

### Initial Guess

- [ ] Gap-aware effective alignment length (P4)

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

## Validation plan

1. **Oscillation test.** Run optimize on datasets known to oscillate (ebola, flu/h3n2/200). Compare likelihood traces with and without damping. Damped runs should show monotonic or near-monotonic likelihood improvement.

2. **Model parity test.** Run optimize with `--model=hky85` on datasets with known transition/transversion bias. Compare branch lengths against ancestral command output with the same model.

3. **Gap initialization test.** Compare initial branch length estimates on gappy alignments (tb, lassa/L) before and after P4. Verify estimates are closer to v0's gap-excluded Hamming distance.

4. **Method equivalence test.** Run optimize with `--opt-method=newton` and `--opt-method=brent` on the same inputs. Final branch lengths should agree within tolerance (1e-4 relative difference). Likelihood should be equal within 1e-6.

5. **Regression test.** Run smoke tests on all datasets with default parameters. Verify no regression in likelihood or branch length output compared to pre-change baseline.

## References

- Brent RP (1973). _Algorithms for Minimization without Derivatives_. Prentice-Hall.
- Clancy D, Lyu S, Roch S (2025). Sample complexity of branch-length estimation by maximum likelihood. arXiv 2507.22038.
- Dinh V, Matsen FA IV (2015). The shape of the one-dimensional phylogenetic likelihood function. arXiv 1507.03647.
- Hasegawa M, Kishino H, Yano T (1985). Dating of the human-ape splitting by a molecular clock of mitochondrial DNA. J Mol Evol 22:160-174.
- Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Mol Biol Evol 37(5):1530-1534.
- Nguyen LT, Schmidt HA, von Haeseler A, Minh BQ (2015). IQ-TREE: A fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Mol Biol Evol 32(1):268-274.
- Polyak BT, Tremba AA (2017). New versions of Newton method: step-size choice, convergence domain and under-determined equations. arXiv 1703.07810.
- Shea D, Schmidt M (2024). Greedy Newton: Newton's method with exact line search. arXiv 2401.06809.
- Stamatakis A (2006). RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models. Bioinformatics 22(21):2688-2690.
- Stamatakis A (2014). RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics 30(9):1312-1313.
