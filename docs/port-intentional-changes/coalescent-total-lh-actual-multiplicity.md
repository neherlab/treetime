# Coalescent total LH uses actual multiplicity instead of fixed 2

v1 uses the parent node's child count for the per-edge multiplicity parameter in the coalescent total log-likelihood formula. v0 uses a fixed default of 2.0 for all edges.

**Type**: Bug fix (v0 erratum correction).

**v0 location**: `MergerModel.total_LH()` at `packages/legacy/treetime/treetime/merger_models.py#L252-L257`. Calls `self.cost(node.time_before_present, node.branch_length)` without passing `multiplicity`, getting the default `multiplicity=2.0`.

**v1 location**: `collect_coalescent_edges()` at `packages/treetime/src/commands/timetree/coalescent/edge_data.rs`. Looks up the parent node via `graph.get_node(parent_node_key)` and uses `parent.outbound().len()` as the multiplicity.

**v0 erratum**: `docs/port-v0-errata/coalescent-total-lh-fixed-multiplicity.md`.

## Background

The per-edge coalescent cost distributes the merger rate credit across the edges entering a parent node:

```
cost = I(t_merger) - I(t_node) - log(λ(t_merger)) * (m-1)/m
```

where m is the number of children at the parent (merger) node. For binary trees (m=2 everywhere), `(m-1)/m = 0.5`. For polytomies, the correct value depends on the parent's child count. v0's `node_contribution()` at `merger_models.py:239-250` correctly uses `multiplicity = len(node.clades)` (parent's child count) for the backward pass inference. `total_LH()` omits the parameter and gets the default 2.0.

## Impact

- Zero impact on fully resolved binary trees (m=2 for all edges)
- For trees with unresolved polytomies: v0 underweights the merger rate credit at polytomy nodes. v1 uses the correct parent-based multiplicity.
- The backward pass inference is unaffected (both v0 and v1 use correct per-node multiplicity there)

## Affected code paths

- `compute_coalescent_total_lh()` - convergence metric
- `TcCostFunction::compute_total_lh()` - Tc optimization objective
- Both call `sum_coalescent_cost()` which uses `CoalescentEdgeData.multiplicity`
