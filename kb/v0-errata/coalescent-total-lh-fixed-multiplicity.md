# total_LH uses fixed multiplicity=2 for all edges

## v0 location

`MergerModel.total_LH()` (`#MergerModel`, `#total_LH`) [packages/legacy/treetime/treetime/merger_models.py#L252-L257](../../packages/legacy/treetime/treetime/merger_models.py#L252-L257)

## Erratum

`total_LH()` calls `self.cost(node.time_before_present, node.branch_length)` without passing the `multiplicity` parameter. `cost()` defaults to `multiplicity=2.0`, treating every edge as a binary merger regardless of actual tree topology. For polytomies (nodes with 3+ children), the merger rate factor `(m-1)/m` is underestimated by 25-33%.

| Children | Correct (m-1)/m | v0 uses (2-1)/2 | Error |
| -------: | --------------: | --------------: | ----: |
|        2 |            0.50 |            0.50 |    0% |
|        3 |            0.67 |            0.50 |  -25% |
|        4 |            0.75 |            0.50 |  -33% |

## Evidence

- `MergerModel.cost()` (`#cost`) at [packages/legacy/treetime/treetime/merger_models.py#L215](../../packages/legacy/treetime/treetime/merger_models.py#L215) accepts a `multiplicity` parameter with docstring: "2 if merger is binary, higher if this is a polytomy"
- `MergerModel.node_contribution()` (`#node_contribution`) at [packages/legacy/treetime/treetime/merger_models.py#L239-L250](../../packages/legacy/treetime/treetime/merger_models.py#L239-L250) correctly uses `multiplicity = len(node.clades)` (actual child count) for the backward pass inference
- `MergerModel.total_LH()` (`#total_LH`) at [packages/legacy/treetime/treetime/merger_models.py#L256](../../packages/legacy/treetime/treetime/merger_models.py#L256) omits the third argument with no comment explaining the choice
- `MergerModel.optimize_Tc()` (`#optimize_Tc`) at [packages/legacy/treetime/treetime/merger_models.py#L259-L279](../../packages/legacy/treetime/treetime/merger_models.py#L259-L279) calls `total_LH()`, so Tc optimization also uses the wrong multiplicity

## v0 impact

- Zero impact on fully resolved binary trees (multiplicity is always 2)
- For trees with unresolved polytomies: `optimize_Tc()` converges to a different value because the objective underweights merger rate terms at polytomy nodes
- Backward pass inference is unaffected (`node_contribution()` uses correct multiplicity)

## v1 status

v1 no longer evaluates the edge-based total likelihood: since phase 2 of
[timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md)
the reported likelihood and the $T_c$ objective are the expected complete-data log-likelihood of the expected
lineage counts, in which each node contributes its actual number of mergers $m_v$
([packages/treetime/src/coalescent/mean_field.rs](../../packages/treetime/src/coalescent/mean_field.rs),
`LineageExpectation::segment_statistics`). The v0 defect therefore has no counterpart in v1.
