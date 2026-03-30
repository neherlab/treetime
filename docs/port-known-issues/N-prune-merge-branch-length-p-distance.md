# Merge shared mutations uses raw p-distance for branch length

`merge_sibling_pair()` in `packages/treetime/src/commands/prune/run.rs` estimates the branch length for the new internal edge as `total_shared_mutations / total_alignment_length`. This is a raw p-distance (Hamming distance), not a model-corrected evolutionary distance.

## Problem

Under JC69, the relationship between observed proportion of differences $p$ and evolutionary distance $d$ is:

$$d = -\frac{3}{4} \ln\left(1 - \frac{4p}{3}\right)$$

The raw ratio underestimates $d$: by 7% at $p = 0.1$, by 18% at $p = 0.25$.

## Impact

Low for the prune command's intended use case. Shared mutation counts between siblings are small relative to alignment length, keeping $p$ small and the approximation error bounded. The child branch length adjustment (`max(0, original - new_edge_bl)`) is also approximate, so correcting only one side has limited benefit.

## Fix

Apply the Jukes-Cantor correction when the model is JC69, with a saturation fallback for $p \geq 0.75$.
