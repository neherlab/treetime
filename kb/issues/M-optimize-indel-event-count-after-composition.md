# Indel event count reduced after composition on merged edges

> [!IMPORTANT]
> **Decision required.** The Poisson indel term is a v1 addition with no v0 counterpart, so v0 cannot settle which count is correct. The options under "Possible resolutions" define different statistics (original event count, net annotation count, or sub-region count) and change the estimated indel rate. Choose one before changing the code.

After `compose_indels()` merges adjacent or overlapping deletions on a collapsed/rerooted edge, `edge_indel_count()` returns the number of net indel annotations (`indels.len()`) rather than the number of indel events that occurred on the original edges.

## Impact

`estimate_indel_rate()` and `total_indel_log_lh()` in `packages/treetime/src/optimize/indel.rs` use `edge_indel_count()` as the Poisson event count `k`. When composition merges two adjacent deletions into one annotation, `k` decreases from 2 to 1 on that edge. This may bias the global indel rate estimate downward.

## Context

The old concatenation approach preserved the raw event count but double-counted positions in `Composition::add_indel()`, producing wrong nucleotide frequency counts for GTR parameter estimation. The composition fix corrects the composition counting (the primary bug) but changes the event count semantics.

## Affected code

- Sparse `fn edge_indel_count()` [packages/treetime/src/partition/marginal/sparse/partition.rs#L99](../../packages/treetime/src/partition/marginal/sparse/partition.rs#L99) returns `indels.len()`; dense `fn edge_indel_count()` [packages/treetime/src/partition/marginal/dense/partition.rs#L154](../../packages/treetime/src/partition/marginal/dense/partition.rs#L154) does the same
- `fn compose_indels()` [packages/treetime/src/seq/indel.rs#L10](../../packages/treetime/src/seq/indel.rs#L10) merges overlapping or adjacent deletions into one annotation
- `fn gather_edge_indel_counts()` [packages/treetime/src/optimize/gather.rs#L55](../../packages/treetime/src/optimize/gather.rs#L55) collects `edge_indel_count()` for every edge
- `fn estimate_indel_rate()` [packages/treetime/src/optimize/indel.rs#L16](../../packages/treetime/src/optimize/indel.rs#L16) sums these counts across all edges
- `fn total_indel_log_lh()` [packages/treetime/src/optimize/indel.rs#L42](../../packages/treetime/src/optimize/indel.rs#L42) uses the count of each edge as the Poisson event count

## Frequency

Only affects edges produced by topology cleanup (collapse) or reroot (merge). Adjacent deletions on consecutive edges require two independent indel events at neighboring positions on adjacent branches, which is uncommon in phylogenetic data.

## Possible resolutions

- Track original event count separately from net annotation count during composition
- Accept annotation count as the Poisson statistic (simpler model, minor bias)
- Weight the Poisson count by the number of indel sub-regions in each annotation
