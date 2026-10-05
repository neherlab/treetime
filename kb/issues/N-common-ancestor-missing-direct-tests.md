# common_ancestor missing direct tests

`fn common_ancestor()` [`packages/treetime-graph/src/common_ancestor.rs#L7-L27`](../../packages/treetime-graph/src/common_ancestor.rs#L7-L27) is only tested indirectly through rerooting callers. Missing direct coverage for: empty input (should error), single-node identity, sibling MRCA, different-depth nodes, and invalid node keys (actionable error, no panic). Tests should assert exact node identities and error classes, not only success or failure.
