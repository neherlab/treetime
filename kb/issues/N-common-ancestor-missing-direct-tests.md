# common_ancestor missing direct tests

`packages/treetime-graph/src/common_ancestor.rs` is only tested indirectly through clock reroot-tips pipeline tests. Missing direct coverage for: empty input (should error), single-node identity, sibling MRCA, different-depth nodes, and invalid node keys.
