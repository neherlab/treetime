# Polytomy resolution discards the edge payload of every reparented child

## Severity

Medium — incorrect behavior under bounded conditions: affects only children reparented during
polytomy resolution, and only when `--resolve-polytomies` is active.

## Location

`fn merge_children()` [packages/treetime/src/timetree/optimization/polytomy.rs#L517-L558](../../packages/treetime/src/timetree/optimization/polytomy.rs#L517-L558)

## Defect

Merging two children under a new internal node reparents each child by removing its edge and
adding a new one with a freshly defaulted payload:

```rust
let new_branch1_length = child1.time - new_node_time;
graph.remove_edge(child1.edge_key)?;
let child1_edge_payload = EdgeTimetree {
  time_length: Some(new_branch1_length),
  ..EdgeTimetree::default()
};
graph.add_edge(new_node_key, child1.node_key, child1_edge_payload)?;
```

Only `time_length` survives. Everything else on `EdgeTimetree`
([packages/treetime/src/payload/timetree.rs#L161-L177](../../packages/treetime/src/payload/timetree.rs#L161-L177))
is reset to its default. Two fields carry state that is still live at this point:

**`branch_length`** — the observed mutation length on the branch, reset to `None`.
`fn prepare_tree_after_topology_change()` runs immediately afterwards and explicitly
documents this field as one of the two it preserves to seed the next inference pass:
"Keep the observed branch length and inferred time length: both seed the next pass"
([packages/treetime/src/timetree/optimization/polytomy.rs#L611-L612](../../packages/treetime/src/timetree/optimization/polytomy.rs#L611-L612)).
That function goes on to null out `branch_length_distribution`, `msg_to_parent`, `gamma` and
the three `ClockSet` fields, which is precisely the set it considers topology-dependent.
`branch_length` is deliberately excluded — and then lost anyway, upstream, for reparented
children only.

**`gamma`** — the relaxed-clock per-branch rate multiplier, reset from its inferred value to
the default `1.0`. `fn apply_relaxed_clock()` assigns it at
[packages/treetime/src/timetree/refinement.rs#L30](../../packages/treetime/src/timetree/refinement.rs#L30),
three lines before topology refinement runs at
[packages/treetime/src/timetree/refinement.rs#L33](../../packages/treetime/src/timetree/refinement.rs#L33).
Under `--relax`, every reparented child silently reverts to the average clock rate while its
unreparented siblings keep theirs.

A third consequence is indirect. `add_edge` allocates a new `GraphEdgeKey`
([packages/treetime-graph/src/graph_ops.rs#L77](../../packages/treetime-graph/src/graph_ops.rs#L77)),
so the partition entry keyed by the old edge is stale. `reconcile_topology` drops it and
inserts a `Default` under the new key
([packages/treetime/src/partition/marginal_sparse.rs#L493-L500](../../packages/treetime/src/partition/marginal_sparse.rs#L493-L500)),
discarding that edge's substitution and message state. The subsequent `update_marginal`
([packages/treetime/src/timetree/refinement.rs#L105](../../packages/treetime/src/timetree/refinement.rs#L105))
recomputes `subs_ml` and the messages, so this one is recovered; `branch_length` and `gamma`
are not.

## Evidence

- `prepare_tree_after_topology_change` enumerates the fields it resets and names
  `branch_length` as intentionally retained, contradicting what `merge_children` does to it
- `gamma` is assigned by the immediately preceding step in the same `Refinement::run()` body
- The parent-side edge created at
  [packages/treetime/src/timetree/optimization/polytomy.rs#L533-L537](../../packages/treetime/src/timetree/optimization/polytomy.rs#L533-L537)
  legitimately wants a default payload — it is a new branch with no history. The child-side
  edges are the same branch as before, relocated, and want their payload carried over.

## Impact

Bounded today: the greedy merger reparents two children per merge and stops once the
likelihood gain falls below `DEFAULT_RESOLUTION_THRESHOLD`, so a typical run touches few
edges. The effect scales with how much of the tree gets reparented, and
[kb/proposals/timetree-stochastic-polytomy-resolution.md](../proposals/timetree-stochastic-polytomy-resolution.md)
would reparent most children of every polytomy.

## Fix

Add `Graph::reparent_edge(edge_key, new_source)` to `treetime-graph`: remove the key from the
old source's `outbound_mut()`, push it to the new source's, and `set_source`. Roughly 25
lines. `Graph::collapse_edge` already reparents this way internally
([packages/treetime-graph/src/graph_ops.rs#L186-L199](../../packages/treetime-graph/src/graph_ops.rs#L186-L199)),
so the mechanism is established; it is simply not exposed.

Keeping the edge key preserves the whole payload and the partition entry, leaving only
`time_length` to overwrite. It also removes an `O(V)` node scan per reparent
(`remove_edge` at
[packages/treetime-graph/src/graph_ops.rs#L106-L124](../../packages/treetime-graph/src/graph_ops.rs#L106-L124)
walks every node).

The same primitive is wanted by
[kb/proposals/optimize-polytomy-reversion-resolution.md](../proposals/optimize-polytomy-reversion-resolution.md),
where it was costed and deferred.

## Test

- Reparent a child whose edge carries a non-default `branch_length` and `gamma`; assert both
  survive on the new edge.
- Under `--relax`, assert that gamma values across a resolved polytomy's children are not all
  `1.0`.
- Assert the edge key is unchanged after reparenting, and that the sparse partition entry for
  that key is preserved.
