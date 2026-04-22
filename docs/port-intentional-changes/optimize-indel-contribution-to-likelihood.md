# Indel contribution to branch length likelihood

## Deviation

v1 includes a Poisson indel term in the per-edge log-likelihood during branch length optimization. v0 ignores indels entirely, treating gaps as missing data.

## Rationale

A branch with zero substitutions but one or more indels represents genuine evolutionary change. Without an indel contribution, such branches are assigned zero length, collapsing topology that the indel evidence supports.

The Poisson model adds $\log P(k \mid \mu t)$ to the edge log-likelihood, where $k$ is the observed indel count, $\mu$ is the global indel rate (estimated as total indels / total branch length), and $t$ is the branch length. The derivatives $k/t - \mu$ and $-k/t^2$ integrate into the existing Newton optimization.

## Impact

Low for most datasets. Indels are rare in typical viral phylogenetics. The effect is visible on branches where the only signal of divergence is an indel event, and on optimize-loop convergence diagnostics for indel-bearing runs. The contribution is a no-op only when all partitions report zero indels.

## Implementation

- `optimize_indel.rs`: Poisson log-likelihood, derivatives, global rate estimation (generic over any `Graph<N, E, ()>` whose edges implement `HasBranchLength`)
- `optimize_unified.rs`: indel contribution added to `run_optimize_mixed`, `initial_guess_mixed`, and the zero-branch optimality check
- `run.rs`: `run_optimize_loop()` records the joint substitution + indel objective and passes one per-iteration `indel_rate` to both tree-level evaluation and per-edge optimization
- `timetree/inference/branch_length_likelihood.rs` and `timetree/inference/runner.rs`: the timetree branch-length distribution grid uses the same `evaluate_with_indels_log_lh_only()` evaluator, with `indel_rate` estimated once per pass and `indel_count` computed per edge
- `partition_ops.rs`: `edge_indel_count()` trait method

## Convergence note

The indel rate $\hat{\mu} = \sum_e k_e / \sum_e t_e$ is estimated from current branch lengths at each optimization round. On the first iteration, branch lengths come from `initial_guess_mixed` which bootstraps indel-only edges to `one_mutation` (a small value). This makes the denominator artificially small and the rate estimate artificially high, biasing branches shorter on the first iteration. The bias self-corrects on subsequent iterations as branch lengths converge.

`run_optimize_loop()` now uses the same per-iteration $\hat{\mu}$ both for the recorded outer-loop likelihood and for `run_optimize_mixed()`. This removes the previous objective mismatch where edge optimization included indels but `LH`, convergence checks, and rollback logic ignored them.

**Update**: investigation of [M-optimize-sparse-em-2-cycle](../port-known-issues/M-optimize-sparse-em-2-cycle.md) confirmed that per-iteration $\hat\mu$ recomputation amplifies a 2-cycle caused by the sparse variable/fixed position boundary. On sc2/2844, $\hat\mu \approx 12{,}000$ (3751 indels / 0.31 total BL), and a 0.06% BL oscillation shifts $\hat\mu$ proportionally across all edges. Proposed fix: compute $\hat\mu$ once before the loop and cache it. See [optimize-convergence-and-robustness](../port-proposals/optimize-convergence-and-robustness.md) P4.

## Double-counting caveat

The indel rate estimator and per-edge count in `run_optimize_mixed()` sum `edge_indel_count()` across all partitions. When dense and sparse partitions represent the same alignment, this produces the correct count only if one partition type has zero indels. Currently, Fitch reconstruction populates indels on sparse partitions only. If indel detection is added for dense partitions, partition-aware deduplication is needed to avoid doubling the count and the Poisson curvature.

## Integration note

`edge_indel_count()` is on `PartitionOptimizeOps`. When `PartitionBranchOps` (branch `refactor/unified-branch-mutations-api`) is merged, move `edge_indel_count` to `PartitionBranchOps` since indel counts are a general partition property.

## Alternatives considered

The Poisson count model was chosen over more sophisticated approaches. See [indel models report](../reports/indel-models/_index.md) for a full catalog of indel modeling approaches with scientific background, and [indel model alternatives proposal](../port-proposals/optimize-indel-model-alternatives.md) for future directions.

Three approaches were evaluated:

1. **Affine gap penalty** - fixed cost per indel event plus per-position extension cost. Not probabilistic; cannot produce a proper likelihood.
2. **TKF91 birth-death process** - separates insertion rate $\lambda$ from deletion rate $\mu$ with equilibrium constraint $\lambda < \mu$. Tracks individual indel lengths. Computationally expensive ($O(L^N)$ exact), requires restructuring the per-edge likelihood. Over-parameterized for the branch-length-prevents-zero use case.
3. **Poisson indel count** (chosen) - single rate, each indel event has equal weight. Negligible computational cost. Integrates directly into Newton step via additive log-likelihood term.

The primary goal is preventing zero-length assignment on branches with only indel evidence, not reconstructing the indel process. The Poisson model achieves this with minimal implementation and computational cost.

## v0 handling

v0 ignores indels in the likelihood. This is consistent with RAxML, IQ-TREE, PhyML, and BEAST, which all treat gaps as missing data. The v1 indel contribution is a design-doc feature, not a v0 port.
