# Indel contribution to branch length likelihood

## Deviation

v1 includes a Poisson indel term in the per-edge log-likelihood during branch length optimization. v0 ignores indels entirely, treating gaps as missing data.

## Rationale

A branch with zero substitutions but one or more indels represents genuine evolutionary change. Without an indel contribution, such branches are assigned zero length, collapsing topology that the indel evidence supports.

The Poisson model adds $\log P(k \mid \mu t)$ to the edge log-likelihood, where $k$ is the observed indel count, $\mu$ is the global indel rate (estimated as total indels / total branch length), and $t$ is the branch length. The derivatives $k/t - \mu$ and $-k/t^2$ integrate into the existing Newton optimization.

## Impact

Low for most datasets. Indels are rare in typical viral phylogenetics. The effect is visible on branches where the only signal of divergence is an indel event, and on optimize-loop convergence diagnostics for indel-bearing runs. The contribution is a no-op only when all partitions report zero indels.

## Implementation

- `optimize/indel.rs`: Poisson log-likelihood, derivatives, global rate estimation
- `optimize/dispatch.rs`: indel contribution added to `run_optimize_mixed`, `initial_guess_mixed`, and the zero-branch optimality check
- `optimize/run_loop.rs`: `run_optimize_loop()` estimates `indel_rate` once from the branch lengths it receives and passes it to both tree-level evaluation and per-edge optimization in every iteration
- `optimize/likelihood.rs`: `evaluate_with_indels_log_lh_only()`, also used by the timetree branch-length distribution grid (`timetree/inference/branch_length_likelihood.rs`, `timetree/inference/runner.rs`), with `indel_rate` estimated once per pass and `indel_count` computed per edge
- `edge_indel_count()` on the sparse and dense partitions and on `MarginalReconstruction`

## Convergence note

The indel rate $\hat{\mu} = \sum_e k_e / \sum_e t_e$ is estimated once at the start of each `run_optimize_loop()` call [`packages/treetime/src/optimize/run_loop.rs#L48-L53`](../../packages/treetime/src/optimize/run_loop.rs#L48-L53) and held fixed for its iterations. It uses the branch lengths that enter the loop, which come from `initial_guess_mixed` when the input lacks valid lengths. That function bootstraps indel-only edges to `one_mutation` (a small value) when no rate is available, which makes the denominator small and the rate estimate high, biasing those branches shorter. Because the rate is fixed during the loop, this bias does not correct itself within one call.

`run_optimize_loop()` uses the same $\hat{\mu}$ both for the recorded outer-loop likelihood and for `run_optimize_mixed()`, so edge optimization, `LH`, convergence checks, and rollback logic all see the same objective.

A rate recomputed in every iteration would amplify a 2-cycle caused by the sparse variable/fixed position boundary. On sc2/2844, $\hat{\mu} \approx 12{,}000$ (3751 indels / 0.31 total BL), and a 0.06% BL oscillation shifts $\hat{\mu}$ proportionally across all edges. Computing $\hat{\mu}$ once before the loop avoids this. See [optimize-convergence-and-robustness](../proposals/optimize-convergence-and-robustness.md) P4.

## Double-counting caveat

The indel rate estimator and per-edge count in `run_optimize_mixed()` sum `edge_indel_count()` across all partitions. When dense and sparse partitions represent the same alignment, this produces the correct count only if one partition type has zero indels. Currently, Fitch reconstruction populates indels on sparse partitions only. If indel detection is added for dense partitions, partition-aware deduplication is needed to avoid doubling the count and the Poisson curvature.

## Integration note

`edge_indel_count()` is on `PartitionOptimizeOps`. Consider moving it to `PartitionBranchOps` since indel counts are a general partition property.

## Alternatives considered

The Poisson count model was chosen over more sophisticated approaches. See the indel models report ([kb/reports/indel-models/1-introduction.md](../reports/indel-models/1-introduction.md)) for a full catalog of indel modeling approaches with scientific background, and [indel model alternatives proposal](../proposals/optimize-indel-model-alternatives.md) for future directions.

Three approaches were evaluated:

1. Affine gap penalty - fixed cost per indel event plus per-position extension cost. Not probabilistic; cannot produce a proper likelihood.
2. TKF91 birth-death process - separates insertion rate $\lambda$ from deletion rate $\mu$ with equilibrium constraint $\lambda < \mu$. Tracks individual indel lengths. Computationally expensive ($O(L^N)$ exact), requires restructuring the per-edge likelihood. Over-parameterized for the branch-length-prevents-zero use case.
3. Poisson indel count (chosen) - single rate, each indel event has equal weight. Negligible computational cost. Integrates directly into Newton step via additive log-likelihood term.

The primary goal is preventing zero-length assignment on branches with only indel evidence, not reconstructing the indel process. The Poisson model achieves this with minimal implementation and computational cost.

## v0 handling

v0 ignores indels in the likelihood. This is consistent with RAxML, IQ-TREE, PhyML, and BEAST, which all treat gaps as missing data. The v1 indel contribution is a v1 addition, not a v0 port.
