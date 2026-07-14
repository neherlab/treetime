# Coalescent leaf contribution comments misstate the node-grouped sign

The coalescent leaf helper returns the correct node-grouped negative-log term, $-I(t_{\mathrm{leaf}})$, but its comments describe that term as an isolated branch-survival probability and then claim that `DistributionNegLog` contains log-probability values. Those descriptions contradict the whole-tree decomposition and obscure the sign contract [packages/treetime/src/coalescent/contributions.rs#L21-L24](../../packages/treetime/src/coalescent/contributions.rs#L21-L24) [packages/treetime/src/coalescent/contributions.rs#L82-L92](../../packages/treetime/src/coalescent/contributions.rs#L82-L92).

## Mathematical contract

Let $I(t)$ be the cumulative per-lineage merger rate at time $t$. For a branch from child $c$ at time $t_c$ to parent $p$ at time $t_p$, the survival contribution to negative log-likelihood is

$$
I(t_p)-I(t_c).
$$

Summing over all branches assigns $-I(t_c)$ to each leaf, $m_j I(t_j)$ to each non-root internal node $j$ with merger multiplicity $m_j$, and $+I(t_{\mathrm{root}})$ once at the root. Therefore the leaf helper's negative value is a node-grouped component of a non-negative whole-tree survival cost; it is neither the negative log-likelihood of an isolated leaf nor evidence that `DistributionNegLog` stores log probabilities. V0 uses the same leaf and root grouping [packages/legacy/treetime/treetime/clock_tree.py#L497-L526](../../packages/legacy/treetime/treetime/clock_tree.py#L497-L526).

## Impact

- The comments invite a sign inversion that would make v1 disagree with both the branch-sum derivation and v0.
- Reviewers cannot infer the `DistributionNegLog` representation contract from the current explanation.
- Existing raw-formula golden masters do not directly state the whole-tree identity that determines the sign.

## Potential solutions

- O1. Correct the node-grouped documentation and add a direct branch-sum identity test.
- O2. Replace node-grouped contributions with explicit per-edge survival terms. This makes every component locally non-negative but changes message-passing ownership and must still reproduce the same whole-tree sum.

## Recommendation

Use O1. Describe the leaf value as a node-grouped negative-log term and explain that it becomes meaningful after internal and root terms are summed. Add a direct unit test comparing the grouped node contributions with the branch-survival sum. Do not change the numerical sign.

## Related issues

- [M-timetree-coalescent-missing-leaf-and-root-contributions.md](M-timetree-coalescent-missing-leaf-and-root-contributions.md)
