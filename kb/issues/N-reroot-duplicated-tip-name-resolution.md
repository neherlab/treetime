# Tip-name resolution duplicated between optimize and clock reroot

> [!IMPORTANT]
> **Discussion required.** The shared lookup contract and duplicate-name semantics need agreement before the two implementations are consolidated.

Two reroot call sites implement the same graph-node-by-name lookup independently. Optimize's `fn resolve_tip_keys()` [packages/treetime/src/optimize/pipeline.rs#L251](../../packages/treetime/src/optimize/pipeline.rs#L251) and clock's `fn find_node_key_by_name()` [packages/treetime/src/clock/reroot.rs#L270](../../packages/treetime/src/clock/reroot.rs#L270) both scan the graph for a node whose name matches a string, using the identical predicate.

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

Both implementations return the first matching storage slot. Duplicate tip names therefore select an arbitrary node while reporting success, and graph construction or topology changes can alter the selected node.

## Potential solutions

- O1. Build one shared name index mapping each name to every matching key and reject ambiguity.
- O2. Require callers to provide stable node keys rather than names. This is precise for library callers but does not solve CLI name input.

## Recommendation

Extract one shared ambiguity-detecting name index and use it from CLI resolution in both reroot paths. The migration of clock rerooting onto the generic reroot module is the natural consolidation point.

## Locations

- `fn resolve_tip_keys()` [packages/treetime/src/optimize/pipeline.rs#L251](../../packages/treetime/src/optimize/pipeline.rs#L251)
- `fn find_node_key_by_name()` [packages/treetime/src/clock/reroot.rs#L270](../../packages/treetime/src/clock/reroot.rs#L270)
- [kb/proposals/reroot-generic-scoring-architecture.md](../proposals/reroot-generic-scoring-architecture.md)
- [kb/issues/N-reroot-tip-resolution-untested-errors.md](N-reroot-tip-resolution-untested-errors.md)
