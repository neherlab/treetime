# Indel contribution to branch length likelihood

## Deviation

v1 includes a Poisson indel term in the per-edge log-likelihood during branch length optimization. v0 ignores indels entirely, treating gaps as missing data.

## Rationale

A branch with zero substitutions but one or more indels represents genuine evolutionary change. Without an indel contribution, such branches are assigned zero length, collapsing topology that the indel evidence supports.

The Poisson model adds $\log P(k \mid \mu t)$ to the edge log-likelihood, where $k$ is the observed indel count, $\mu$ is the global indel rate (estimated as total indels / total branch length), and $t$ is the branch length. The derivatives $k/t - \mu$ and $-k/t^2$ integrate into the existing Newton optimization.

## Impact

Low for most datasets. Indels are rare in typical viral phylogenetics. The effect is visible on branches where the only signal of divergence is an indel event. For standard datasets (flu, ebola, zika), the indel rate is zero and the contribution is a no-op.

## Implementation

- `optimize_indel.rs`: Poisson log-likelihood, derivatives, global rate estimation
- `optimize_unified.rs`: indel contribution added to `run_optimize_mixed`, `initial_guess_mixed`, and the zero-branch optimality check
- `partition_ops.rs`: `edge_indel_count()` trait method

## v0 handling

v0 ignores indels in the likelihood. This is consistent with RAxML, IQ-TREE, PhyML, and BEAST, which all treat gaps as missing data. The v1 indel contribution is a design-doc feature, not a v0 port.
