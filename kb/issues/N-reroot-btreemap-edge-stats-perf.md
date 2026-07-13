# Reroot edge-stats maps use BTreeMap where Vec indexing fits

`packages/treetime/src/reroot/div_stats_traversal.rs` stores directional DivStats messages in four `BTreeMap<GraphEdgeKey, _>` maps. `GraphEdgeKey(pub usize)` is a dense index -- the graph stores edges in `Vec<Option<Edge>>`. Replacing with `Vec<Option<DivStats>>` indexed by `edge_key.as_usize()` would give O(1) access instead of O(log E). Once-per-optimize-run, low priority. Same applies to `find_best_split` at `packages/treetime/src/reroot/split.rs:66`.
