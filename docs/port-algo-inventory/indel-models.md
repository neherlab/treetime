# Indel Models

[Back to index](_index.md)

For the full survey of indel modeling approaches in phylogenetics (12 models, 17 software tools, 30 references), see the [Indel Models in Phylogenetics](../reports/indel-models/_index.md) report.

This page covers only the algorithm implemented in v1.

## Poisson Indel Contribution

Adds a Poisson indel log-likelihood term to per-edge branch length optimization. For $k$ observed indel events on a branch of length $t$ with global rate $\mu$: $\ell(t) = k \ln(\mu t) - \mu t - \ln(k!)$. Derivatives $k/t - \mu$ and $-k/t^2$ enter the Newton step alongside substitution derivatives. The rate $\hat{\mu} = \sum_e k_e / \sum_e t_e$ is estimated from the tree at each optimization round.

v1: [`packages/treetime/src/commands/optimize/optimize_indel.rs`](../../packages/treetime/src/commands/optimize/optimize_indel.rs).

v0: not implemented. v0 ignores indels in the likelihood, same as RAxML, IQ-TREE, PhyML.

This is a v1-only feature. See [indel models report](../reports/indel-models/_index.md) for the full catalog of alternative approaches, [intentional change](../port-intentional-changes/optimize-indel-contribution-to-likelihood.md), and [alternatives proposal](../port-proposals/optimize-indel-model-alternatives.md).
