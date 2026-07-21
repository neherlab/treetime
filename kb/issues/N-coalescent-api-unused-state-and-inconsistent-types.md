# Coalescent APIs retain inconsistent types and unused return state

The coalescent module's public API has type inconsistencies and unused state that complicate maintenance and hide potential bugs.

## Unused present-time return

`fn collect_tree_events()` [`packages/treetime/src/coalescent/events.rs#L14`](../../packages/treetime/src/coalescent/events.rs#L14) returns `(CalendarTime, Vec<...>)` where the `CalendarTime` is the present time (maximum leaf date). The sole caller [`packages/treetime/src/coalescent/precomputed.rs#L36`](../../packages/treetime/src/coalescent/precomputed.rs#L36) discards it with `let (_, events)`.

## Inconsistent child-count representations

`struct CoalescentEdge` stores `n_children: f64` [`packages/treetime/src/coalescent/edge_data.rs#L16`](../../packages/treetime/src/coalescent/edge_data.rs#L16), while `fn CoalescentModel::internal_contribution()` and `fn root_contribution()` accept `n_children: usize` [`packages/treetime/src/coalescent/coalescent.rs#L65-L72`](../../packages/treetime/src/coalescent/coalescent.rs#L65-L72). The edge collector converts `outbound().len()` (usize) to `f64` via `as` at [`packages/treetime/src/coalescent/edge_data.rs#L99-L100`](../../packages/treetime/src/coalescent/edge_data.rs#L99-L100), while the model methods convert the other direction with `as f64` inside arithmetic. Neither direction is checked.

## Repeated parent lookup with unreachable fallback

Already documented in [M-coalescent-edge-collection-nan-bypass-and-unreachable-fallback.md](M-coalescent-edge-collection-nan-bypass-and-unreachable-fallback.md) (section "Repeated parent lookup with unreachable fallback"). Included here for completeness of the API-level view.

## Required behavior

Use a single typed representation for child counts (`usize` at construction, converted once to `f64` where arithmetic requires it). Remove unused return values. Remove redundant lookups and unreachable fallbacks.

## Related issues

- [M-coalescent-edge-collection-nan-bypass-and-unreachable-fallback.md](M-coalescent-edge-collection-nan-bypass-and-unreachable-fallback.md)
