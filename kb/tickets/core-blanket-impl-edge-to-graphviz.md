# Add blanket or default implementation for EdgeToGraphviz trait

Four edge types implement EdgeToGraphviz with identical bodies that format branch_length.

## Current state

`EdgeTimetree`, `EdgeAncestral`, `EdgeClock`, and `TestEdge` each implement `to_graphviz_label()` and `to_graphviz_weight()` identically. All implement `HasBranchLength`.

## Target state

Either a blanket `impl<T: HasBranchLength> EdgeToGraphviz for T` or default methods on the `EdgeToGraphviz` trait itself. Types needing custom graphviz rendering override.

## Implementation

1. Check if `EdgeToGraphviz` is defined in a crate that can see `HasBranchLength` (for blanket impl) or if default methods are more appropriate
2. Add default method implementations using `self.branch_length()` and `format_weight`
3. Remove the four identical explicit implementations
4. Verify graphviz output unchanged for all edge types

## Related issues

Source: `kb/issues/N-core-edge-to-graphviz-identical-impls.md`
