# Edge annotations from Newick not wired into the parsed graph

`graph_from_newick()` in `packages/treetime-io/src/nwk.rs` keeps only the branch length of each edge (`NwkParse.branch_lengths`). Parsed branch-level comments and annotations (`NewickEdgeData::comments()`) from `util-newick` are discarded.

BEAST2 canonical format places branch annotations after `:` (e.g. `A:[&rate=0.003]0.1`). These are parsed by `util-newick` into the edge comments but never reach `NwkParse`.

No current command reads or writes branch-level annotations. All written annotation data comes from `fn nwk_node_comments()` in `packages/app-output/src/nwk_comments.rs`, which returns node-keyed `NwkNodeComments` (`packages/treetime-io/src/nwk.rs`). Wiring would require an edge-keyed annotation map in `NwkParse` beside `branch_lengths`, and an edge-keyed comment map for the writers.

## Locations

- Discarded data: `packages/treetime-io/src/nwk.rs` `graph_from_newick()` edge loop
- Parse result: `packages/treetime-io/src/nwk.rs` `NwkParse`
- Parsed data: `packages/util-newick/src/model/data.rs` `NewickEdgeData::comments()`

## Related issues

- The writer emits annotations in node position (before `:`), not branch position (after `:`), so round-trip through treetime preserves node annotations without needing edge annotation wiring

## Proposed implementation plan

1. Parse: add `pub branch_annotations: EdgeMap<Vec<(String, NewickValue)>>` to `struct NwkParse` (`packages/treetime-io/src/nwk.rs`). In the edge loop of `fn graph_from_newick()`, insert `nwk_edge.data().annotations()` (`NewickEdgeData::annotations()` in `packages/util-newick/src/model/data.rs`, which joins the occurrence comments and the edge comments) beside `branch_lengths`
2. Consumer first: wire the map only when a command needs it (for example BEAST2 `rate` input). Without a reader the field is dead code that `hawk` reports, so this step and step 1 land together
3. Write (only if a command must write branch annotations): add an edge-keyed `NwkEdgeComments` beside `NwkNodeComments` and emit it after `:` in the Newick and Nexus writers (`packages/treetime-io/src/nwk.rs`, `packages/treetime-io/src/nex.rs`), filled from `packages/app-output/src/nwk_comments.rs`
4. Tests in `packages/treetime-io/src/__tests__/test_nwk_comments.rs`: parse `(A:[&rate=0.003]0.1,B:0.2);` and compare the whole edge annotation map with the expected literal (oracle: BEAST2 Newick annotation syntax, value from the input); round-trip test when step 3 is done
5. Validate: `just lint-rs`, `just test-rs -E 'package(treetime-io)'`
