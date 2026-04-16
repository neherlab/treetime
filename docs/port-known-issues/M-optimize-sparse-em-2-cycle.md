# Optimize outer loop enters persistent likelihood 2-cycle (sparse)

## Problem

The sparse optimize loop fails to converge on datasets with indels. The log-likelihood alternates between two values (a "2-cycle") on stable topology, with amplitude exceeding the convergence threshold `dp`. The loop exhausts all `max_iter` iterations.

Dense mode converges correctly on the same dataset.

## Reproduction

```bash
./dev/docker/run ./dev/dev r treetime -- optimize \
  --tree data/sc2/2844/tree.nwk --outdir tmp/oscillation_test \
  --dense false --aln data/sc2/2844/aln.fasta.xz \
  --branch-length-initial-guess always --model jc69 --max-iter 50 -vv
```

Expected: convergence before `max_iter`. Actual: likelihood alternates ~-143156/-143157 from iteration ~15 onward. Tree topology is stable from iteration 5.

## Root cause

The sparse representation classifies alignment positions as "variable" (individual posterior distribution per position) or "fixed" (canonical-state profile shared by all same-state positions, counted by multiplicity). The sparse evaluator `get_coefficients()` (`#get_coefficients`) [packages/treetime/src/commands/optimize/optimize_sparse.rs#L45](../../packages/treetime/src/commands/optimize/optimize_sparse.rs#L45) computes eigendecomposition coefficients from individual posteriors for variable positions and from bucketed canonical profiles for fixed positions. These two paths produce different coefficients for the same position because the individual posterior incorporates the full tree context while the canonical profile does not.

The sparse marginal forward/backward passes (`process_node_backward()`, `process_node_forward()`) (`#process_node_backward`, `#process_node_forward`) [packages/treetime/src/representation/partition/marginal_passes.rs#L34](../../packages/treetime/src/representation/partition/marginal_passes.rs#L34) rebuild the variable position set from scratch each iteration. On sc2/2844, ~14 positions oscillate between variable and fixed classification, confirmed by per-iteration instrumentation:

```
Iteration 11: total_var_msg_entries=35095   LH=-1.431562e5
Iteration 12: total_var_msg_entries=35109   LH=-1.431568e5
Iteration 13: total_var_msg_entries=35096   LH=-1.431563e5
Iteration 14: total_var_msg_entries=35108   LH=-1.431559e5
```

The variable count alternates in lockstep with the likelihood. On datasets where the variable set is stable (flu/h3n2/20: 1121 every iteration), the loop converges normally.

Dense mode [packages/treetime/src/representation/partition/marginal_dense.rs](../../packages/treetime/src/representation/partition/marginal_dense.rs) has no variable/fixed classification. Every position stores a full probability vector in an `Array2<f64>` matrix. The dense evaluator [packages/treetime/src/commands/optimize/optimize_dense.rs](../../packages/treetime/src/commands/optimize/optimize_dense.rs) computes coefficients from these matrices directly. The E-step is purely continuous. The dense `reconstruct_node_sequence()` (`#reconstruct_node_sequence`) [packages/treetime/src/representation/partition/marginal_dense.rs#L341](../../packages/treetime/src/representation/partition/marginal_dense.rs#L341) uses `argmax_first()` (`#argmax_first`) for output sequences, but `update_marginal()` (`#update_marginal`) [packages/treetime/src/commands/ancestral/marginal.rs#L41](../../packages/treetime/src/commands/ancestral/marginal.rs#L41) -- the function called during iteration -- invokes only `marginal_backward()` + `marginal_forward()` (`#marginal_backward`, `#marginal_forward`) [packages/treetime/src/commands/ancestral/marginal.rs#L85-L122](../../packages/treetime/src/commands/ancestral/marginal.rs#L85-L122), not `reconstruct_node_sequence()`. Dense mode is true soft-EM with guaranteed monotone likelihood increase <a id="cite-1"></a>([Dempster, Laird, and Rubin 1977](https://doi.org/10.1111/j.2517-6161.1977.tb01600.x) [[1](#ref-1)]).

The sparse variable/fixed boundary is an inherent approximation of the sparse representation. The discrete reclassification produces a non-monotone objective, violating the monotone convergence guarantee of standard EM <a id="cite-4"></a>([Wu 1983](https://doi.org/10.1214/aos/1176346060) [[4](#ref-4)]). The convergence loop must be robust to this non-monotonicity rather than assuming the smooth EM convergence guarantee. The sparse evaluator effectively operates as a generalized ECM <a id="cite-2"></a>([Meng and Rubin 1993](https://doi.org/10.1093/biomet/80.2.267) [[2](#ref-2)]) with a discrete conditional maximization step (the variable/fixed reclassification), but the conditional step is non-smooth, so the ECM monotonicity guarantee does not apply either.

## Amplifying factors

Two factors worsen the 2-cycle but do not cause it:

**Indel rate feedback loop.** The global indel rate $\hat{\mu} = \sum_e k_e / \sum_e t_e$ is recomputed each iteration in `run_optimize_mixed()` (`#run_optimize_mixed`) [packages/treetime/src/commands/optimize/optimize_unified.rs#L532](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L532), where $k_e$ is the observed indel count on edge $e$ and $t_e$ is the branch length. On sc2/2844 with 3751 indels and ~0.31 total BL, $\hat{\mu} \approx 12{,}000$. A 0.06% BL oscillation shifts $\hat{\mu}$ proportionally, amplifying the discrete jump across all edges via the Poisson term $k \ln(\hat{\mu} t) - \hat{\mu} t - \ln(k!)$.

**Damping decay.** Exponential damping $d^{i+1}$ (default $d = 0.75$) decays geometrically: ~0.042 at iteration 10, ~0.003 at iteration 20. By the time the 2-cycle emerges (~iteration 15), damping is insufficient to bridge the discrete jump.

## Contributing code issues

**Commented-out `estimate_indel_rate()` (`#estimate_indel_rate`) in `initial_guess_mixed()` (`#initial_guess_mixed`).** PR [neherlab/treetime#558](https://github.com/neherlab/treetime/pull/558) replaced the call with `let indel_rate = 0.0;` at [packages/treetime/src/commands/optimize/optimize_unified.rs#L718](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L718). Three experimental runs on sc2/2844 show restoring the call is the only variant that converges:

| Run              | `initial_guess_mixed`        | `run_optimize_mixed`           | Converges?    | Final LH        | 2-cycle amplitude |
| ---------------- | ---------------------------- | ------------------------------ | ------------- | --------------- | ----------------- |
| A (current code) | hardcoded `0.0`              | `estimate_indel_rate` (~12000) | NO (50 iter)  | -143156/-143157 | ~1.4              |
| B (restored)     | `estimate_indel_rate` (~9.3) | `estimate_indel_rate` (~12000) | YES (iter 17) | -143151         | 0                 |
| C (both `0.0`)   | hardcoded `0.0`              | hardcoded `0.0`                | NO (30 iter)  | -142143/-142148 | ~0.05             |

Run B converges because overestimated initial branch lengths force slow descent through damped iterations. By the time the optimizer nears the optimum, damping is still ~1.3% ($0.75^{15}$), bridging the discrete jump.

**Unjustified v1 default mismatch.** v0 uses `max_iter=10` and `LHtol=0.1` in `optimize_tree_marginal()` (`#optimize_tree_marginal`) [packages/legacy/treetime/treetime/treeanc.py#L1297-L1298](../../packages/legacy/treetime/treetime/treeanc.py#L1297-L1298). v1 uses `max_iter=20` and `dp=0.01` [packages/treetime/src/commands/optimize/args.rs#L129-L134](../../packages/treetime/src/commands/optimize/args.rs#L129-L134). No documented reason for the divergence. v1's tighter threshold and higher iteration count push the loop into the undamped regime.

Note: v0's convergence check is signed (`deltaLH < LHtol`), accepting any likelihood decrease as convergence. v1 uses absolute value. The signed check masks non-convergence as an accidental side effect. See [optimize-signed-convergence-check](../port-v0-errata/optimize-signed-convergence-check.md) for the erratum.

## Fix

Four changes that make the convergence loop robust to the sparse non-monotone objective. No changes to the sparse marginal passes or the optimization model.

### 1. Restore `estimate_indel_rate` in `initial_guess_mixed`

Revert the PR #558 change: replace `let indel_rate = 0.0;` with `let indel_rate = estimate_indel_rate(graph, partitions);`.

Location: `initial_guess_mixed()` (`#initial_guess_mixed`) [packages/treetime/src/commands/optimize/optimize_unified.rs#L700](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L700).

### 2. Three-condition convergence check

Replace the single abs-delta check with three orthogonal stopping conditions in `run_optimize_loop()` (`#run_optimize_loop`) [packages/treetime/src/commands/optimize/run.rs#L260](../../packages/treetime/src/commands/optimize/run.rs#L260), where $\text{LH}_i$ is the log-likelihood at iteration $i$ and $\mathit{dp}$ is the convergence threshold:

- **Converged**: $|\text{LH}_i - \text{LH}_{i-1}| < \mathit{dp}$ (current criterion, catches true convergence)
- **Oscillating**: $|\text{LH}_i - \text{LH}_{i-2}| < \mathit{dp}$ (detects 2-cycles where consecutive differences exceed $\mathit{dp}$ but the two-step difference is small). Requires $i \ge 2$.
- **Worsened**: $\text{LH}_i < \text{LH}_{i-1}$ and $\text{LH}_{i-1}$ was the best observed likelihood so far (likelihood decreased after a peak; restore branch lengths from the best-observed state and stop). Requires $i \ge 2$.

The worsened condition with revert ensures the output always reflects the best observed state. This is similar to IQ-TREE's per-round monotonicity check with rollback <a id="cite-3"></a>([Minh et al. 2020](https://doi.org/10.1093/molbev/msaa015) [[3](#ref-3)]).

Implementation: ~20 lines. Maintain `best_lh` and `best_branch_lengths` updated whenever a new best likelihood is observed. Store `lh_prev_prev` alongside `lh_prev` for oscillation detection. The oscillation and worsened checks are skipped for the first 2 iterations (insufficient history). On the worsened condition, restore `best_branch_lengths` before returning.

### 3. Damping floor

Replace `pow(damping, iteration + 1)` with `pow(damping, iteration + 1).max(DAMPING_FLOOR)` where `DAMPING_FLOOR = 0.01` in `apply_damping()` (`#apply_damping`) [packages/treetime/src/commands/optimize/run.rs#L366](../../packages/treetime/src/commands/optimize/run.rs#L366). Prevents fully undamped late iterations regardless of `max_iter`.

### 4. Match v0 defaults

Change `max_iter` default from 20 to 10 and `dp` default from 0.01 to 0.1 in [packages/treetime/src/commands/optimize/args.rs#L129-L134](../../packages/treetime/src/commands/optimize/args.rs#L129-L134). Matching v0's `optimize_tree_marginal()` (`#optimize_tree_marginal`) defaults keeps the loop in the well-damped regime for typical usage.

### Cleanup

Remove investigation diagnostic `eprintln!` blocks from [packages/treetime/src/commands/optimize/optimize_unified.rs](../../packages/treetime/src/commands/optimize/optimize_unified.rs) (`[DIAG run_optimize_mixed]` and `[DIAG initial_guess_mixed]`), if still present.

## Verification

### Reproduce the bug

Build current `dev` and run the reproduction command. Verify: likelihood 2-cycle from iteration ~15, no convergence.

### Verify the fix

1. Apply all four changes. Run on sc2/2844 with `--max-iter 50`. Verify: loop converges (via any of the three conditions).
2. Run with `--max-iter 100`. Verify: no late-onset oscillation.
3. Run on flu/h3n2/20 (no indels, stable variable set). Verify: convergence unchanged.
4. Run in dense mode on flu/h3n2/20. Verify: convergence unchanged.

### Regression tests

- **Convergence on sc2/2844**: integration test running `run_optimize_loop` (sparse, JC69, `--branch-length-initial-guess always`, max_iter=50), assert `stopped_at.is_some()`.
- **Damping floor**: unit test verifying `apply_damping` retains nonzero old-value weight at iteration 1000.
- **Oscillation detection**: unit test with synthetic likelihood sequence exhibiting a 2-cycle, verify the oscillation condition triggers.
- **Revert-on-decrease**: unit test verifying branch lengths are restored to previous iteration when likelihood decreases.
- **Convergence on indel-free dataset**: flu/h3n2/20 converges within 10 iterations.

## v1 status

**Fixed.** All four changes implemented:

1. `estimate_indel_rate` restored in `initial_guess_mixed()` [packages/treetime/src/commands/optimize/optimize_unified.rs#L718](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L718)
2. Three-condition convergence check (`ConvergenceReason::Converged`, `Oscillating`, `Worsened`) in `run_optimize_loop()` [packages/treetime/src/commands/optimize/run.rs#L260](../../packages/treetime/src/commands/optimize/run.rs#L260)
3. `DAMPING_FLOOR = 0.01` in `apply_damping()` [packages/treetime/src/commands/optimize/run.rs#L383](../../packages/treetime/src/commands/optimize/run.rs#L383)
4. Defaults aligned with v0: `max_iter=10`, `dp=0.1` [packages/treetime/src/commands/optimize/args.rs#L129-L134](../../packages/treetime/src/commands/optimize/args.rs#L129-L134)

Diagnostic `eprintln!` blocks were already absent at implementation time.

Verified on sc2/2844 (the original reproduction dataset, sparse JC69, converges within 50 iterations), flu/h3n2/20 (sparse JC69, converges within 10 iterations), and the toy tree test suite (10 unit tests covering all three conditions, damping floor, and branch length revert).

## Related

- [M-optimize-gm-per-branch-divergence](M-optimize-gm-per-branch-divergence.md) -- per-branch v0 parity; non-convergence may contribute
- [M-optimize-iterative-log-likelihood-nan](M-optimize-iterative-log-likelihood-nan.md) -- NaN on larger datasets; sparse forward-pass instability
- [M-gtr-sparse-composition-stale-after-marginal](M-gtr-sparse-composition-stale-after-marginal.md) -- stale Fitch composition in GTR inference
- [optimize-signed-convergence-check](../port-v0-errata/optimize-signed-convergence-check.md) -- v0 erratum: signed convergence check conflates convergence with failure
- [optimize-convergence-and-robustness](../port-proposals/optimize-convergence-and-robustness.md) -- proposal for deeper architectural improvements (variable-set stabilization, evaluator fixes, convergence acceleration)
- [docs/algorithms/optimize.md](../algorithms/optimize.md) -- design document
- [docs/reports/iterative-tree-refinement/9-iteration-loop.md](../reports/iterative-tree-refinement/9-iteration-loop.md) -- damping scheme, convergence theory

## References

1. <a id="ref-1"></a> Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. "Maximum Likelihood from Incomplete Data via the EM Algorithm." _Journal of the Royal Statistical Society: Series B_ 39 (1): 1-38. [doi:10.1111/j.2517-6161.1977.tb01600.x](https://doi.org/10.1111/j.2517-6161.1977.tb01600.x) [↩](#cite-1)
2. <a id="ref-2"></a> Meng, Xiao-Li, and Donald B. Rubin. 1993. "Maximum Likelihood Estimation via the ECM Algorithm: A General Framework." _Biometrika_ 80 (2): 267-278. [doi:10.1093/biomet/80.2.267](https://doi.org/10.1093/biomet/80.2.267) [↩](#cite-2)
3. <a id="ref-3"></a> Minh, B. Q., H. A. Schmidt, O. Chernomor, D. Schrempf, M. D. Woodhams, A. von Haeseler, and R. Lanfear. 2020. "IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era." _Molecular Biology and Evolution_ 37 (5): 1530-1534. [doi:10.1093/molbev/msaa015](https://doi.org/10.1093/molbev/msaa015) [↩](#cite-3)
4. <a id="ref-4"></a> Wu, C. F. Jeff. 1983. "On the Convergence Properties of the EM Algorithm." _The Annals of Statistics_ 11 (1): 95-103. [doi:10.1214/aos/1176346060](https://doi.org/10.1214/aos/1176346060) [↩](#cite-4)
