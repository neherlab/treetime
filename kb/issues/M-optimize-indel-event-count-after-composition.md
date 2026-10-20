# Indel event count reduced after composition on merged edges

> [!IMPORTANT]
> **Decision required.** The Poisson indel term is a v1 addition with no v0 counterpart, so v0 cannot settle which count is correct. The options under "Possible resolutions" define different statistics (original event count, net annotation count, or sub-region count) and change the estimated indel rate. Choose one before changing the code.

In sparse partitions, after `compose_indels()` merges adjacent or overlapping deletions of two edges into one edge, `edge_indel_count()` returns the number of net indel annotations (`indels.len()`) rather than the number of indel events that occurred on the original edges. Dense partitions are not affected: their count comes from the indels of the current marginal estimate of each edge, which no edge merge composes.

## Impact

`estimate_indel_rate()` and `total_indel_log_lh()` in `packages/treetime/src/optimize/indel.rs` use `edge_indel_count()` as the Poisson event count `k`. When composition merges two adjacent deletions into one annotation, `k` decreases from 2 to 1 on that edge. This may bias the global indel rate estimate downward.

## Context

Edge merges compose the indel lists instead of concatenating them, because a concatenated list counts overlapping positions twice in `Composition::add_indel()` and gives wrong nucleotide frequency counts for GTR parameter estimation. Concatenation would keep the raw event count; composition keeps the frequency counts correct but changes the event count semantics.

## Affected code

- Sparse `fn edge_indel_count()` [packages/treetime/src/partition/marginal/sparse/partition.rs#L103-L105](../../packages/treetime/src/partition/marginal/sparse/partition.rs#L103-L105) returns `indels.len()` of the stored edge annotations. Dense `fn edge_indel_count()` [packages/treetime/src/partition/marginal/dense/partition.rs#L157-L163](../../packages/treetime/src/partition/marginal/dense/partition.rs#L157-L163) returns `indels.len()` of the marginal estimate of the edge
- `fn compose_indels()` [packages/treetime/src/seq/indel.rs#L9](../../packages/treetime/src/seq/indel.rs#L9) merges overlapping or adjacent deletions into one annotation. It returns the other list unchanged when one list is empty
- `fn gather_edge_indel_counts()` [packages/treetime/src/optimize/gather.rs#L87](../../packages/treetime/src/optimize/gather.rs#L87) collects `edge_indel_count()` for every edge
- `fn estimate_indel_rate()` [packages/treetime/src/optimize/indel.rs#L155](../../packages/treetime/src/optimize/indel.rs#L155) sums these counts across all edges
- `fn total_indel_log_lh()` [packages/treetime/src/optimize/indel.rs#L181](../../packages/treetime/src/optimize/indel.rs#L181) uses the count of each edge as the Poisson event count

## Frequency

Composition merges annotations only when both edges carry indels. These paths compose two such edges:

- Reroot: `SparseEdgeObs::chain_fitch_indels()` on the merged parent and child edge [packages/treetime/src/partition/marginal/sparse/reroot.rs#L112](../../packages/treetime/src/partition/marginal/sparse/reroot.rs#L112)
- Reversion hoist during polytomy resolution: `fn split_indels()` composes the parent and child indels when they interact [packages/treetime/src/optimize/topology/hoist_reversions.rs#L300-L319](../../packages/treetime/src/optimize/topology/hoist_reversions.rs#L300-L319)

The collapse of zero-length edges in the optimize loop and the retirement of polytomy helper edges collapse only edges without indels, so `compose_indels()` returns the child list unchanged there [packages/treetime/src/optimize/run_loop.rs#L312-L331](../../packages/treetime/src/optimize/run_loop.rs#L312-L331), [packages/treetime/src/optimize/topology/resolve_polytomy.rs#L163-L169](../../packages/treetime/src/optimize/topology/resolve_polytomy.rs#L163-L169). Adjacent deletions on consecutive edges require two independent indel events at neighboring positions on adjacent branches, which is uncommon in phylogenetic data.

## Possible resolutions

- Track original event count separately from net annotation count during composition
- Accept annotation count as the Poisson statistic (simpler model, minor bias)
- Weight the Poisson count by the number of indel sub-regions in each annotation
