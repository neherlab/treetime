# Indel event count reduced after composition on merged edges

After `compose_indels()` merges adjacent or overlapping deletions on a collapsed/rerooted edge, `edge_indel_count()` returns the number of net indel annotations (`indels.len()`) rather than the number of indel events that occurred on the original edges.

## Impact

`estimate_indel_rate()` and `total_indel_log_lh()` in `packages/treetime/src/optimize/indel.rs` use `edge_indel_count()` as the Poisson event count `k`. When composition merges two adjacent deletions into one annotation, `k` decreases from 2 to 1 on that edge. This may bias the global indel rate estimate downward.

## Context

The old concatenation approach preserved the raw event count but double-counted positions in `Composition::add_indel()`, producing wrong nucleotide frequency counts for GTR parameter estimation. The composition fix corrects the composition counting (the primary bug) but changes the event count semantics.

## Affected code

- `fn edge_indel_count()` at `packages/treetime/src/partition/marginal_sparse.rs` returns `indels.len()`
- `fn estimate_indel_rate()` at `packages/treetime/src/optimize/indel.rs` sums `edge_indel_count()` across all edges
- `fn total_indel_log_lh()` at `packages/treetime/src/optimize/indel.rs` uses `edge_indel_count()` per edge

## Frequency

Only affects edges produced by topology cleanup (collapse) or reroot (merge). Adjacent deletions on consecutive edges require two independent indel events at neighboring positions on adjacent branches, which is uncommon in phylogenetic data.

## Possible resolutions

- Track original event count separately from net annotation count during composition
- Accept annotation count as the Poisson statistic (simpler model, minor bias)
- Weight the Poisson count by the number of indel sub-regions in each annotation
