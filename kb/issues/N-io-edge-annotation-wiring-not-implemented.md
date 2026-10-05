# Edge annotations from Newick not wired into the parsed graph

`graph_from_newick()` in `packages/treetime-io/src/nwk.rs` keeps only the branch length of each edge (`NwkParse.branch_lengths`). Parsed branch-level annotations (`NewickEdgeData.branch_attrs`) from `util-newick` are discarded.

BEAST2 canonical format places branch annotations after `:` (e.g. `A:[&rate=0.003]0.1`). These are parsed by `util-newick` into `branch_attrs` but never reach `NwkParse`.

No current command reads or writes branch-level annotations. All written annotation data comes from `fn nwk_node_comments()` in `packages/app-output/src/nwk_comments.rs`, which returns node-keyed `NwkNodeComments` (`packages/treetime-io/src/nwk.rs`). Wiring would require an edge-keyed annotation map in `NwkParse` beside `branch_lengths`, and an edge-keyed comment map for the writers.

## Locations

- Discarded data: `packages/treetime-io/src/nwk.rs` `graph_from_newick()` edge loop
- Parse result: `packages/treetime-io/src/nwk.rs` `NwkParse`
- Parsed data: `packages/util-newick/src/types.rs` `NewickEdgeData.branch_attrs`

## Related issues

- The writer emits annotations in node position (before `:`), not branch position (after `:`), so round-trip through treetime preserves node annotations without needing edge annotation wiring
