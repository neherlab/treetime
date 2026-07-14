# Tree-backed output order is inconsistent

Commands apply topology ordering separately from TreeIR construction. Direct formats consume the ordered graph, while most TreeIR-backed formats project the original graph, so the same requested output plan can contain different node orders.

## Evidence

Affected command runners apply topology order under [`packages/treetime/src/commands`](../../packages/treetime/src/commands), while TreeIR projections are constructed separately under [`packages/treetime-io/src/tree_ir`](../../packages/treetime-io/src/tree_ir). `ancestral`, `clock`, `mugration`, `optimize`, and `prune` can therefore ignore the requested topology order for Auspice, PhyloXML, or UShER while direct Newick/Nexus output follows it.

## Potential solutions

- O1. Build one ordered graph and derive every direct output and TreeIR projection from it.
- O2. Reapply the selected order independently inside every writer. This duplicates ordering and can drift as writers are added.

## Recommendation

Use O1. Apply topology order once during output preparation and pass the resulting graph to every writer or projection. Method/format availability is an independent contract tracked by [N-ancestral-auspice-json-not-produced.md](N-ancestral-auspice-json-not-produced.md); this issue does not select a parsimony output policy.

## Related issues

- [N-io-tree-ir-architecture-unapproved.md](N-io-tree-ir-architecture-unapproved.md)
- [N-ancestral-auspice-json-not-produced.md](N-ancestral-auspice-json-not-produced.md)
