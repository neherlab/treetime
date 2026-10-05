# Auspice entropy output perturbs the Shannon definition

> [!IMPORTANT]
> **Decision required.** The perturbed formula is exact Augur parity. Augur computes the trait entropy as `S = -np.sum(pdis*np.log(pdis+TINY))` with `TINY = 1e-12` (`augur/traits.py` lines 13 and 87 at the pinned revision `d8e38736037ba9474a809f9a5a63bc2b279d2407`). Augur is a reference for TreeTime output, so exact Shannon entropy diverges from it and needs approval and a `kb/decisions/` entry. The options are:
>
> - Keep Augur parity and close this issue
> - Use exact Shannon entropy as described under Fix, and record the divergence from Augur

Entropy is computed with an added `TINY` inside every logarithm. This changes every positive term and gives a deterministic distribution such as $(1,0)$ a nonzero entropy.

`fn compute_entropy()` applies $p\ln(p+10^{-12})$ to every state [packages/app-output/src/mugration_tree_output.rs#L177-L180](../../packages/app-output/src/mugration_tree_output.rs#L177-L180). Both the mugration node-data JSON [packages/app-output/src/augur_node_data_mugration.rs#L121](../../packages/app-output/src/augur_node_data_mugration.rs#L121) and the Auspice tree-output entropy [packages/app-output/src/mugration_tree_output.rs#L126](../../packages/app-output/src/mugration_tree_output.rs#L126) consume this function.

For state probabilities $p_i$, Shannon entropy is

$$H=-\sum_i p_i\log p_i$$

where $H$ is entropy, $p_i$ is the probability of state $i$, and a zero-probability term contributes zero by continuity. The logarithm is natural, matching the existing output contract.

## Potential solutions

- O1. Use `ndarray-stats::EntropyExt` with explicit invalid-input errors. `ndarray-stats` is a workspace dependency of `treetime`; `app-output` must add it.
- O2. Implement the $p_i=0$ limit directly in project code. This duplicates a maintained dependency already in use.

## Fix

If exact Shannon entropy is approved:

- Take one borrowed profile row as input
- Use `EntropyExt::entropy()` without adding an epsilon to the probabilities
- Return errors for empty, non-finite, and invalid profiles, with node and trait context

## Validation

- Deterministic, uniform, zero-containing, empty, and invalid distributions, with independent analytical expected values
- Projection test over the whole Auspice document

## Related issues

- [M-timetree-tree-output-inference-metadata-incomplete.md](M-timetree-tree-output-inference-metadata-incomplete.md)
- [N-mugration-confidence-rows-copied-for-output.md](N-mugration-confidence-rows-copied-for-output.md): changes the `compute_entropy()` input to a borrowed row view
