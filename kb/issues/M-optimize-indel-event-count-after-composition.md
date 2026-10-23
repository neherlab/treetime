# Indel event count reduced after composition on merged edges

> [!IMPORTANT]
> **Decision required.** The Poisson indel term is a v1 addition with no v0 counterpart, so v0 cannot settle which count is correct. The options under "Possible resolutions" define different statistics (original event count, net annotation count, or sub-region count) and change the estimated indel rate. Choose one before changing the code.

In sparse partitions, after `compose_indels()` merges adjacent or overlapping deletions of two edges into one edge, the indel count of the merged edge is the number of net indel annotations of its alignment rather than the number of indel events that occurred on the original edges. Composition runs per alignment on the assembled indels of both edges (`fn AlignmentLayout::compose_indels()`), so a split alignment gives the same count as the unsplit alignment. Dense partitions are not affected: their count comes from the indels of the current marginal estimate of each edge, which no edge merge composes.

## Impact

`estimate_indel_rate()` and `total_indel_log_lh()` in `packages/treetime/src/optimize/indel.rs` use the per-alignment indel count of each edge as the Poisson event count `k`. When composition merges two adjacent deletions into one annotation, `k` decreases from 2 to 1 on that edge. This may bias the global indel rate estimate downward.

## Context

Edge merges compose the indel lists instead of concatenating them, because a concatenated list counts overlapping positions twice in `Composition::add_indel()` and gives wrong nucleotide frequency counts for GTR parameter estimation. Concatenation would keep the raw event count; composition keeps the frequency counts correct but changes the event count semantics.

## Affected code

- Sparse `fn edge_indel_count()` [packages/treetime/src/partition/marginal/sparse/partition.rs#L104-L106](../../packages/treetime/src/partition/marginal/sparse/partition.rs#L104-L106) returns `indels.len()` of the stored edge annotations. Dense `fn edge_indel_count()` [packages/treetime/src/partition/marginal/dense/partition.rs#L158-L164](../../packages/treetime/src/partition/marginal/dense/partition.rs#L158-L164) returns `indels.len()` of the marginal estimate of the edge
- `fn compose_indels()` [packages/treetime/src/seq/indel.rs#L9](../../packages/treetime/src/seq/indel.rs#L9) merges overlapping or adjacent deletions into one annotation. It returns the other list unchanged when one list is empty
- `fn gather_alignment_edge_indel_counts()` [packages/treetime/src/optimize/gather.rs#L98-L120](../../packages/treetime/src/optimize/gather.rs#L98-L120) counts the indel events of every edge per alignment: `edge_indel_count()` for a whole alignment, the length of the assembled list (`fn assemble_indels()`) for a split alignment. `fn gather_edge_indel_counts()` [packages/treetime/src/optimize/gather.rs#L88-L96](../../packages/treetime/src/optimize/gather.rs#L88-L96) sums these counts over the alignments
- `fn estimate_indel_rate()` [packages/treetime/src/optimize/indel.rs#L155](../../packages/treetime/src/optimize/indel.rs#L155) sums these counts across all edges
- `fn total_indel_log_lh()` [packages/treetime/src/optimize/indel.rs#L181](../../packages/treetime/src/optimize/indel.rs#L181) uses the count of each edge as the Poisson event count

## Frequency

Composition merges annotations only when both edges carry indels. These paths compose two such edges:

- Edge merges: `fn merge_edge_pair()` [packages/treetime/src/partition/marginal/sparse/alignment_indels.rs#L76-L95](../../packages/treetime/src/partition/marginal/sparse/alignment_indels.rs#L76-L95) composes the parent and child edge per alignment with `fn composed_edge_pieces()` [packages/treetime/src/partition/marginal/sparse/alignment_indels.rs#L26-L36](../../packages/treetime/src/partition/marginal/sparse/alignment_indels.rs#L26-L36). Two paths call it: the merge of the trivial old root in a sparse reroot (`fn remove_trivial_root()` in `packages/treetime/src/partition/marginal/sparse/reroot.rs`), and the edge merges of `fn PartitionSet::reconcile_topology()` through `fn merge_sparse_edges()`
- Reversion hoist during polytomy resolution: `fn split_indels()` composes the assembled parent and child indels of each alignment when they interact [packages/treetime/src/optimize/topology/hoist_reversions.rs#L310-L329](../../packages/treetime/src/optimize/topology/hoist_reversions.rs#L310-L329)

`fn collapse_edge()` also composes per alignment with `composed_edge_pieces()`, but the collapse of zero-length edges in the optimize loop and the retirement of polytomy helper edges collapse only edges without substitutions and indels, so `compose_indels()` returns the child list unchanged there [packages/treetime/src/optimize/run_loop.rs#L334-L353](../../packages/treetime/src/optimize/run_loop.rs#L334-L353), [packages/treetime/src/optimize/topology/resolve_polytomy.rs#L184-L189](../../packages/treetime/src/optimize/topology/resolve_polytomy.rs#L184-L189). Adjacent deletions on consecutive edges require two independent indel events at neighboring positions on adjacent branches, which is uncommon in phylogenetic data.

## Possible resolutions

- Track original event count separately from net annotation count during composition
- Accept annotation count as the Poisson statistic (simpler model, minor bias)
- Weight the Poisson count by the number of indel sub-regions in each annotation
