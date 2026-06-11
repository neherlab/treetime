# Edge annotations from Newick not wired into EdgeFromNwk

`EdgeFromNwk::from_nwk(weight: Option<f64>)` receives only the branch length. Parsed branch-level annotations (`NewickEdgeData.branch_attrs`) from `util-newick` are discarded during `graph_from_newick()` at `packages/treetime-io/src/nwk.rs`.

BEAST2 canonical format places branch annotations after `:` (e.g. `A:[&rate=0.003]0.1`). These are parsed by `util-newick` into `branch_attrs` but never passed to `treetime-io` edge payloads.

No current command reads or writes branch-level annotations -- all annotation data flows through `NodeToNwk::nwk_comments()` and `CommentProviders`, which are node-scoped. Wiring would require extending the `EdgeFromNwk` trait signature to accept a `BTreeMap<String, String>` parameter, updating all implementations.

## Locations

- Trait: `packages/treetime-io/src/nwk.rs` `EdgeFromNwk::from_nwk`
- Discarded data: `packages/treetime-io/src/nwk.rs` `graph_from_newick()` edge loop
- Parsed data: `packages/util-newick/src/types.rs` `NewickEdgeData.branch_attrs`

## Related issues

- Resolved by ticket #1: writer now emits annotations in node position (before `:`), not branch position (after `:`), so round-trip through treetime preserves node annotations without needing edge annotation wiring
