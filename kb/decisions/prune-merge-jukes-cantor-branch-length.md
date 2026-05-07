# Merged-sibling branch length uses Jukes-Cantor correction, not raw p-distance

The `prune --merge-shared-mutations` step groups sibling branches in a polytomy that share substitutions under a new internal node and assigns the new edge a branch length. v1 assigns the Jukes-Cantor 1969 corrected evolutionary distance rather than the raw Hamming ratio specified in the v1 design document.

## What the spec says

[do../_raw/optimize.md](../_raw/optimize.md) proposes the feature for `optimize`:

> scanning for shared subs in children of a polytomy and introducing a new internal node with the individuals that share a sub as children. Initial branch length for this new internal node would be `#shared subs/length`.

`#shared subs / length` is the raw p-distance. Under any reversible substitution model this underestimates the expected number of substitutions per site, because repeated or back-mutations at the same site mask earlier changes.

## What v1 does

[`merge_sibling_pair()`](../../packages/treetime/src/commands/prune/run.rs#L483) applies the Jukes-Cantor 1969 correction through [`jukes_cantor_distance()`](../../packages/treetime/src/gtr/jc_distance.rs#L50):

$$ d = -\frac{k-1}{k} \ln\!\left(1 - \frac{k}{k-1}\, p\right) $$

where $p$ is the pooled p-distance (shared substitutions divided by total alignment length across all partitions) and $k$ is the alphabet size read from the first partition. The formula matches the JC69 model currently hardcoded in [`run_prune()`](../../packages/treetime/src/commands/prune/run.rs#L48) and generalises naturally to amino-acid alphabets (k=20) without code changes.

Saturation at $p \to (k-1)/k$ is handled by clamping $p$ at $p_{sat}\,(1 - 10^{-6})$ before applying the formula. This keeps the result finite (about 10 substitutions per site for k=4, 13 for k=20), continuous in $p$, and safe from `log(0)` or NaN downstream.

## Why v1 differs

Raw p-distance systematically underestimates evolutionary distance. Under JC69 the bias reaches 7% at $p = 0.1$ and 22% at $p = 0.25$. The corrected value is the maximum-likelihood estimate of the branch length under the model the prune command is already assuming (JC69), so applying it at edge creation time is strictly better than deferring the correction to a later ML optimization pass.

The child-branch adjustment `max(0, original - new_edge_bl)` still clamps at zero, so the sum of the new internal edge and a child edge may differ slightly from the child's original length when the correction pushes the internal edge longer than the child. This is unchanged from the raw p-distance behaviour: the same clamp applied before the fix.

## Affected commands

- [`prune --merge-shared-mutations`](../../packages/treetime/src/commands/prune/run.rs#L108) - direct user-facing effect
- [`optimize`](../../packages/treetime/src/commands/optimize/run.rs#L487) - calls `merge_shared_mutation_branches` in its topology-cleanup step before every per-edge optimization pass

## Tests

- Unit: [`test_jukes_cantor_distance_*`](../../packages/treetime/src/gtr/jc_distance.rs#L87) cover known analytical values, monotonicity, the $d \ge p$ property, small-$p$ Taylor behaviour, and the saturation clamp.
- Integration: [`test_merge_branch_length_jc_correction_differs_from_raw`](../../packages/treetime/src/commands/prune/__tests__/test_merge_shared_mutations.rs#L578) exercises a polytomy where $p = 0.1$, regressing if raw p-distance were restored.
