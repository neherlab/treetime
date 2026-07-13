# Shared TreeIR architecture has not received an explicit design decision

> [!IMPORTANT]
> **Research and discussion required.** TreeIR is implemented, but its status as an accepted v1 architecture is unresolved. No decision entry may be added without explicit approval.

## Problem

TreeIR provides one format-neutral graph for PhyloXML, Auspice v2, and UShER MAT serialization. Commands project domain graphs and analysis results into this representation before dispatching format-specific writers. This is a new v1 architecture rather than direct parity with a v0 abstraction.

The implementation defines a shared loss and ownership boundary: domain-specific data must be projected into `TreeIrNode`, `TreeIrEdge`, and `TreeIrData`, while formats outside TreeIR serialize directly from the domain graph. [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md) describes a common graph intermediate representation, but it predates the concrete TreeIR boundary, refers to the former `ConverterGraph`, and does not define TreeIR's preservation or loss invariants. The decision also requires a consent audit under the current project rule that architecture decisions must have explicit approval.

## Research required

- Compare the data requirements and loss behavior of each TreeIR-backed format.
- Establish invariants for node identity, topology order, branch data, mutations, dates, traits, and round trips.
- Decide whether TreeIR is an internal implementation detail or a project-level architectural commitment.
- If the architecture is accepted, obtain explicit approval before adding a `kb/decisions/` entry. If it is not accepted, specify the replacement boundary before implementation work.

## Locations

- `type TreeIrGraph` [packages/treetime-io/src/graph.rs#L23](../../packages/treetime-io/src/graph.rs#L23)
- `pub fn write_tree_outputs()` [packages/treetime-io/src/graph.rs#L61](../../packages/treetime-io/src/graph.rs#L61)
- TreeIR adapters [packages/treetime-io/src/tree_ir/mod.rs](../../packages/treetime-io/src/tree_ir/mod.rs)
- Command projection [packages/treetime/src/commands/shared/ir_projection.rs](../../packages/treetime/src/commands/shared/ir_projection.rs)

## Related KB items

- [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md)
- [kb/proposals/output-format-selection.md](../proposals/output-format-selection.md)
- [kb/proposals/tree-format-model-embedding.md](../proposals/tree-format-model-embedding.md)
- [kb/issues/N-io-auspice-number-formatting-missing-augur-golden-master.md](N-io-auspice-number-formatting-missing-augur-golden-master.md)
