# Extract coalescent event setup into CoalescentPrecomputed struct

Four coalescent functions duplicate 5 lines for tree event collection and lineage count computation.

## Current state

`coalescent.rs:70-76`, `total_lh.rs:38-44`, `optimize_tc.rs:104-110`, `skyline.rs:86-92` each call `collect_tree_events` -> time conversion -> `compute_lineage_count_distribution`.

## Target state

A `CoalescentPrecomputed` struct bundles `present_time`, `events_tbp`, and `lineage_counts` with a `from_graph(graph) -> Result<Self>` constructor. All four sites create the struct and destructure or borrow fields.

## Implementation

1. Define `CoalescentPrecomputed { present_time, events_tbp, lineage_counts }` in `coalescent/` (e.g., `coalescent/precomputed.rs` or in an existing module)
2. Implement `CoalescentPrecomputed::from_graph(graph) -> Result<Self>` with the 5-line setup
3. Replace the 4 call sites with `CoalescentPrecomputed::from_graph(graph)?`
4. Each site accesses fields via the struct (some need `lineage_counts`, some also need `events_tbp`)

## Related issues

Source: `kb/issues/N-coalescent-event-setup-duplication.md`
