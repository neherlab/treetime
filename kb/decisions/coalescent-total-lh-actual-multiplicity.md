# Coalescent total LH uses actual multiplicity instead of fixed 2

v1 uses the parent node's child count for the per-edge multiplicity parameter in the coalescent total log-likelihood formula. v0 uses a fixed default of 2.0 for all edges.

> [!WARNING]
> **Needs review.** The edge-based total likelihood this decision describes (`collect_coalescent_edges`,
> `compute_coalescent_total_lh`, `CoalescentModel`) was removed in phase 2 of
> [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md).
> The reported coalescent likelihood and the convergence metric are now the expected complete-data log-likelihood
> $-\sum_i I_i/T_{c,i} - \sum_i M_i \ln T_{c,i}$, where each node contributes its actual number of mergers $m_v$.

**Type**: Bug fix (v0 erratum correction).

**v0 location**: `MergerModel.total_LH()` at `packages/legacy/treetime/treetime/merger_models.py#L252-L257`. Calls `self.cost(node.time_before_present, node.branch_length)` without passing `multiplicity`, getting the default `multiplicity=2.0`.

**v1 location**: `LineageExpectation::from_survivals()` at `packages/treetime/src/coalescent/mean_field.rs#L341-L357` gives each node the weight $1 - m_v$ from its child count `node.outbound().len()` (L354), so a node with $m_v$ children counts $m_v - 1$ mergers. `LineageExpectation::segment_statistics()` at `packages/treetime/src/coalescent/mean_field.rs#L486-L520` sums these merger counts per skyline segment.

**v0 erratum**: `../v0-errata/coalescent-total-lh-fixed-multiplicity.md`.

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

- `CoalescentTimescale::expected_log_likelihood()` at `packages/treetime/src/timetree/coalescent_timescale.rs#L115` - convergence metric and reported coalescent likelihood
- `optimize_skyline()` at `packages/treetime/src/coalescent/skyline.rs#L14` - Tc optimization objective
- Both use the merger counts of `LineageExpectation::segment_statistics()` in `expected_log_likelihood()` at `packages/treetime/src/coalescent/skyline.rs#L93`
- `CoalescentPrior::node_contribution()` at `packages/treetime/src/coalescent/mean_field.rs#L114` takes the child count `n_children` for the per-node factor of the backward pass
