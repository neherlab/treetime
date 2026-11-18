# common_ancestor missing direct tests

`fn common_ancestor()` [`packages/treetime-graph/src/common_ancestor.rs#L7-L27`](../../packages/treetime-graph/src/common_ancestor.rs#L7-L27) is only tested indirectly through rerooting callers. Missing direct coverage for: empty input (should error), single-node identity, sibling MRCA, different-depth nodes, and invalid node keys (actionable error, no panic). Tests should assert exact node identities and error classes, not only success or failure.

The function has no test in `packages/treetime-graph/src/__tests__/`. Its only callers are `packages/treetime/src/reroot/named_node.rs` and the midpoint test helper in `packages/treetime/src/reroot/__tests__/test_midpoint.rs`. An unknown key returns the error of `Graph::path_from_root_to_node()` (`packages/treetime-graph/src/graph.rs`), which is an internal report ("Node not found on the graph: ..."), not a user-facing error.

## Proposed implementation plan

1. Add `packages/treetime-graph/src/__tests__/test_common_ancestor.rs` and declare it in `packages/treetime-graph/src/__tests__/mod.rs`. Build the graph with the `named_tree` helper pattern of `test_clades.rs` (nodes `root`, `AB`, `A`, `B`, `C`, `D` below `C` for a deeper tip), so the tests name nodes, not keys
2. One `rstest` table, cases named by input class, expected node name from the hand-drawn topology (oracle: the definition of the most recent common ancestor on that tree):
   - single node `[A]` -> `A`
   - siblings `[A, B]` -> `AB`
   - different depths `[A, D]` -> `root`
   - ancestor and descendant `[AB, A]` -> `AB`
   - duplicate keys `[A, A]` -> `A`
   - root included `[root, B]` -> `root`
3. Error cases with `assert_error!` and the exact message: empty slice -> "Cannot find MRCA of an empty node set"; a key that is not in the graph -> the message of `path_from_root_to_node()`
4. Decide whether an unknown key is a caller bug (keep `make_internal_report!`) or a user input error (named-node rerooting passes keys that come from names, so it can not produce an unknown key). Keep the internal report unless a caller can pass user-derived keys, and record the reason in the commit message
5. Validate: `just test-rs -E 'package(treetime-graph)'`, then `just lint-rs`. See red by changing one expected node name and confirming the failure
