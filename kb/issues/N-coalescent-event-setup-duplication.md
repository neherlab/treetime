# Coalescent event setup repeated at four sites

Four coalescent functions repeat a 5-line setup: `collect_tree_events(graph)` -> calendar-to-TBP time conversion -> `compute_lineage_count_distribution`.

## Locations

- `packages/treetime/src/coalescent/coalescent.rs:70-76:` `compute_coalescent_contributions()`
- `packages/treetime/src/coalescent/total_lh.rs:38-44:` `compute_coalescent_total_lh()`
- `packages/treetime/src/coalescent/optimize_tc.rs:104-110:` `TcCostFunction::new()`
- `packages/treetime/src/coalescent/skyline.rs:86-92:` `optimize_skyline()`

Identical pattern at all sites:

```rust
let (present_time, events_calendar) = collect_tree_events(graph)?;
let events_tbp: Vec<_> = events_calendar
  .iter()
  .map(|(t, delta)| (t.to_tbp(present_time), *delta))
  .collect();
let lineage_counts = compute_lineage_count_distribution(&events_tbp)?;
```

## Impact

Pure maintainability. If the event collection or time conversion logic changes, four sites need updating.

## Action

Extract a `CoalescentPrecomputed` struct with a `from_graph(graph) -> Result<Self>` constructor that bundles `present_time`, `events_tbp`, and `lineage_counts`.
