# Coalescent edge collection bypasses NaN dates and retains unreachable multiplicity fallback

Three validation defects in `collect_coalescent_edges`:

## NaN endpoint dates bypass the ordering guard

At [`packages/treetime/src/coalescent/edge_data.rs#L85-L89`](../../packages/treetime/src/coalescent/edge_data.rs#L85-L89), the guard `if child_time < parent_time` uses `f64` comparison. When either time is NaN, the comparison returns `false`, so the edge passes through without error. NaN dates then propagate into the coalescent model's `expected_mergers.eval()` and `total_merger_rate()`, producing silent NaN in the log-likelihood.

## Repeated parent lookup with unreachable fallback

At [`packages/treetime/src/coalescent/edge_data.rs#L92-L93`](../../packages/treetime/src/coalescent/edge_data.rs#L92-L93):

```rust
let n_children = graph
    .get_node(parent_node_key)
    .map_or(2.0, |parent| parent.read_arc().outbound().len() as f64);
```

The parent was already successfully retrieved at lines 69-83. The second `get_node` call is redundant, and the `map_or(2.0, ...)` fallback to binary multiplicity silently masks a graph invariant violation that should be impossible. The fallback is unreachable in correct execution but would produce wrong results if triggered.

## Tc validation deficiency

`CoalescentModel::total_merger_rate()` validates Tc at each queried coordinate ([`coalescent.rs#L84-L98`](../../packages/treetime/src/coalescent/coalescent.rs#L84-L98)), but construction-time validation occurs only at lineage-count midpoints during integration ([`integration.rs#L42-L59`](../../packages/treetime/src/coalescent/integration.rs#L42-L59)). A Distribution with valid midpoints but invalid endpoint or interior values would pass construction and fail only at later query time, with a less informative error context.

## Related issues

- [N-timetree-negative-coalescent-tc.md](N-timetree-negative-coalescent-tc.md): negative Tc accepted without validation (related but distinct -- that issue covers CLI-level acceptance, not model construction deficiencies)
- [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md): topology violations are the upstream source of some reversed dates that trigger the ordering guard
