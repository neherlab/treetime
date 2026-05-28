# EdgeToGraphviz trait has four identical implementations

Four edge types implement `EdgeToGraphviz` with identical method bodies: `to_graphviz_label()` formats `branch_length()` via `format_weight`, `to_graphviz_weight()` returns `branch_length()`. All types already implement `HasBranchLength`.

## Locations

- `packages/treetime/src/payload/timetree.rs:262:` `impl EdgeToGraphviz for EdgeTimetree`
- `packages/treetime/src/payload/ancestral.rs:94:` `impl EdgeToGraphviz for EdgeAncestral`
- `packages/treetime/src/clock/clock_graph.rs:106:` `impl EdgeToGraphviz for EdgeClock`
- `packages/treetime/src/test_utils/graph_fixtures.rs:67:` `impl EdgeToGraphviz for TestEdge`

## Impact

Pure maintainability. New edge types require copying the same boilerplate.

## Action

Add a blanket `impl<T: HasBranchLength> EdgeToGraphviz for T` or provide default method implementations on the trait itself. Types needing custom graphviz rendering can override.
