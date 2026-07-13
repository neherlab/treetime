# Tip-name resolution duplicated between optimize and clock reroot

Two reroot call sites implement the same graph-node-by-name lookup independently. Optimize's `resolve_tip_keys` and clock's `find_node_key_by_name` both scan the graph for a node whose name matches a string, using the identical predicate.

## Duplication

Optimize (`packages/treetime/src/optimize/pipeline.rs:251-260`):

```rust
fn resolve_tip_keys(graph: &GraphAncestral, tips: &[String]) -> Result<Vec<GraphNodeKey>, Report> {
  tips.iter().map(|tip| {
    graph.find_node(|node| node.name().is_some_and(|name| name.as_ref() == tip))
      .ok_or_else(|| make_report!("Reroot tip not found: {tip}"))
  }).collect()
}
```

Clock (`packages/treetime/src/clock/reroot.rs:274-281`):

```rust
fn find_node_key_by_name<N, E, D>(graph: &Graph<N, E, D>, name: &str) -> Option<GraphNodeKey>
where N: GraphNode + Named, E: GraphEdge, D: Send + Sync {
  graph.find_node(|node| node.name().is_some_and(|node_name| node_name.as_ref() == name))
}
```

The core lookup is identical; the optimize version adds multi-tip mapping and a not-found error.

## Impact

Negligible. Two call sites, stable behavior. The duplication is a maintenance hazard: a future change to name-matching semantics (case handling, whitespace, disambiguation of duplicate names) must be applied in both places.

## Fix

Extract a single name-to-node-key lookup (generic over `GraphNode + Named`), place it where both reroot paths can share it (the graph crate or a shared reroot helper), and have `resolve_tip_keys` build on it. The migration of clock rerooting onto the generic reroot module (`kb/tickets/reroot-migrate-clock-to-generic-search.md`) is a natural point to consolidate.

## Locations

- `packages/treetime/src/optimize/pipeline.rs:251-260`
- `packages/treetime/src/clock/reroot.rs:274-281`
