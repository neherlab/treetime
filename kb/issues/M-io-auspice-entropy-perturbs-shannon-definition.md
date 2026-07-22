# Auspice entropy output perturbs the Shannon definition

Entropy is computed with an added `TINY` inside every logarithm. This changes every positive term and gives a deterministic distribution such as $(1,0)$ a nonzero entropy.

`fn compute_entropy()` applies $p\ln(p+10^{-12})$ to every state [packages/treetime/src/commands/mugration/augur_node_data.rs#L135-L137](../../packages/treetime/src/commands/mugration/augur_node_data.rs#L135-L137). Both the mugration node-data JSON [packages/treetime/src/commands/mugration/augur_node_data.rs#L113](../../packages/treetime/src/commands/mugration/augur_node_data.rs#L113) and the Auspice tree-output entropy [packages/treetime/src/commands/shared/tree_output.rs#L1385](../../packages/treetime/src/commands/shared/tree_output.rs#L1385) consume this function.

For state probabilities $p_i$, Shannon entropy is

$$H=-\sum_i p_i\log p_i$$

where $H$ is entropy, $p_i$ is the probability of state $i$, and a zero-probability term contributes zero by continuity. The logarithm is natural, matching the existing output contract.

## Potential solutions

- O1. Use `ndarray-stats::EntropyExt` with explicit invalid-input errors.
- O2. Implement the $p_i=0$ limit directly in project code. This duplicates a maintained dependency already in use.

## Recommendation

Use `ndarray-stats::EntropyExt`, already present in the dependency graph, and propagate empty or invalid input errors. Do not perturb valid probabilities to avoid defining the zero term.

## Related issues

- [M-timetree-tree-output-inference-metadata-incomplete.md](M-timetree-tree-output-inference-metadata-incomplete.md)
