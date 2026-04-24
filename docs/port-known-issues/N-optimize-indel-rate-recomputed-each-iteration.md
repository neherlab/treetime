# Indel rate is recomputed each optimization iteration

## Problem

`estimate_indel_rate()` is called inside `compute_iteration_likelihood()` on every iteration of the optimize loop. The global indel rate $\hat{\mu} = \sum_e k_e / \sum_e t_e$ is recomputed from current branch lengths and passed to both the likelihood evaluation and the per-edge branch-length optimizer.

The GTR substitution rate `gtr.mu` is a comparable global rate parameter but is estimated once before the loop and held fixed. The indel rate does not follow the same convention.

## Motivation

Computing the indel rate once before the loop would:

- Match the treatment of `gtr.mu` (estimated once, held fixed)
- Remove a feedback path where branch-length changes shift $\hat\mu$, which shifts branch-length targets, which shifts $\hat\mu$ again
- Reduce amplification of the sparse 2-cycle on large datasets (documented in P4 of the [convergence proposal](../port-proposals/optimize-convergence-and-robustness.md))

On sc2/2844, $\hat\mu \approx 12{,}000$ (3751 indels / 0.31 total branch length) and oscillates proportionally with total branch length across iterations.

## Affected code

| Location                                                                                                                                   | Role                                                                   |
| ------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------- |
| [packages/treetime/src/commands/optimize/optimize_indel.rs#L55](../../packages/treetime/src/commands/optimize/optimize_indel.rs#L55)       | `estimate_indel_rate()` definition                                     |
| [packages/treetime/src/commands/optimize/run.rs#L397](../../packages/treetime/src/commands/optimize/run.rs#L397)                           | Call site inside `compute_iteration_likelihood()`                      |
| [packages/treetime/src/commands/optimize/run.rs#L359](../../packages/treetime/src/commands/optimize/run.rs#L359)                           | `run_optimize_mixed_with_indel_rate()` consumes the per-iteration rate |
| [packages/treetime/src/commands/optimize/run.rs#L389](../../packages/treetime/src/commands/optimize/run.rs#L389)                           | `compute_iteration_likelihood()` recomputes rate each call             |
| [packages/treetime/src/commands/optimize/optimize_unified.rs#L552](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L552) | `run_optimize_mixed_with_indel_rate()` definition                      |

## Solutions

### S1. Hoist before the loop

Compute `estimate_indel_rate()` once before the first iteration. Pass the cached value to both `compute_iteration_likelihood()` and `run_optimize_mixed_with_indel_rate()`.

- Minimal change (~10 lines)
- Removes the feedback path
- Changes the optimization objective (fixed-rate model vs coordinate-ascent model)

### S2. Cache after topology stabilizes

Recompute while `prune_and_merge_in_loop()` returns `true` (topology still changing, indel counts may change). Cache once topology stabilizes.

- Preserves correct $\hat\mu$ during topology changes when indel counts shift
- Same steady-state behavior as S1
- Proposed as P4 in the [convergence proposal](../port-proposals/optimize-convergence-and-robustness.md)

### S3. Keep current behavior, fix the root cause

Address the sparse variable/fixed boundary discontinuity (P1-P3 in the [convergence proposal](../port-proposals/optimize-convergence-and-robustness.md)). The 2-cycle that $\hat\mu$ recomputation amplifies is caused by the boundary, not by the rate update itself.

## Disagreement

There is a counterargument that per-iteration re-estimation is the theoretically correct ECM approach and that `gtr.mu` being fixed is itself a known missing feature (P9 in the convergence proposal), not a precedent to extend. See the full analysis in the [indel rate re-estimation report](../reports/optimize-indel-rate-reestimation.md).
