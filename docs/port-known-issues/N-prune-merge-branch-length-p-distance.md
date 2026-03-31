# Merge shared mutations uses raw p-distance for branch length

`merge_sibling_pair()` in `packages/treetime/src/commands/prune/run.rs` estimates the branch length for the new internal edge as `total_shared_mutations / total_alignment_length`. This is a raw p-distance (Hamming distance), not a model-corrected evolutionary distance.

## Problem

Under JC69, the relationship between observed proportion of differences $p$ and evolutionary distance $d$ is <a id="cite-1"></a>[Jukes and Cantor 1969](https://doi.org/10.1016/B978-1-4832-3211-9.50009-7) [[1](#ref-1)]:

$$d = -\frac{3}{4} \ln\left(1 - \frac{4p}{3}\right)$$

The raw ratio underestimates $d$: by 7% at $p = 0.1$, by 18% at $p = 0.25$.

## v0 comparison

v0 does not have shared-mutation merging (`merge_shared_mutation_branches` is v1-only). The design spec in [docs/algorithms/optimize.md](../algorithms/optimize.md) says "Initial branch length for this new internal node would be #shared subs/length", which is the raw p-distance currently implemented. The JC69 correction is an improvement over the spec.

## Impact

Low for the prune command's intended use case. Shared mutation counts between siblings are small relative to alignment length, keeping $p$ small and the approximation error bounded. The child branch length adjustment (`max(0, original - new_edge_bl)`) is also approximate, so correcting only one side has limited benefit.

## Related

- [docs/algorithms/optimize.md](../algorithms/optimize.md) -- design spec defining the merge-shared-mutations feature

## Fix

Apply the Jukes-Cantor correction when the model is JC69, with a saturation fallback for $p \geq 0.75$.

## References

1. <a id="ref-1"></a> Jukes, Thomas H., and Charles R. Cantor. 1969. "Evolution of Protein Molecules." In _Mammalian Protein Metabolism_, vol. 3, edited by Hans N. Munro, 21-132. Academic Press. https://doi.org/10.1016/B978-1-4832-3211-9.50009-7 [↩](#cite-1)
