# Shared TreeIR architecture has not received an explicit design decision

> [!IMPORTANT]
> **Research and discussion required.** TreeIR is implemented, but its status as an accepted v1 architecture is unresolved. No decision entry may be added without explicit approval.

## Problem

TreeIR provides one format-neutral graph for PhyloXML, Auspice v2, and UShER MAT serialization. Commands project domain graphs and analysis results into this representation before dispatching format-specific writers. This is a new v1 architecture rather than direct parity with a v0 abstraction.

The implementation defines a shared loss and ownership boundary: domain-specific data must be projected into `TreeIrNode`, `TreeIrEdge`, and `TreeIrData`, while formats outside TreeIR serialize directly from the domain graph. [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md) describes a common graph intermediate representation, but it predates the concrete TreeIR boundary, refers to the former `ConverterGraph`, and does not define TreeIR's preservation or loss invariants. The decision also requires a consent audit under the current project rule that architecture decisions must have explicit approval.

## Decision axes

### A1. Architectural boundary

- O1. Keep TreeIR as an internal, format-neutral projection boundary. Commands produce one typed intermediate value and writers consume only the fields they declare. This centralizes topology and analysis-data projection without making TreeIR a public compatibility promise.
- O2. Project each domain graph directly into each output format. Each writer can model its schema precisely, but ordering, mutation conversion, and typed-field validation would be repeated across adapters.
- O3. Adopt TreeIR as a project-level architecture with a stable preservation contract. This gives every producer and consumer one durable target, but every new format requirement must remain representable in that shared contract.

**Recommendation:** O1. Retain one internal projection boundary and document its invariants, while allowing its types to change when a format exposes a missing domain concept. This preserves centralized conversion without freezing an unshipped internal API.

### A2. Preservation and unsupported states

- O1. Model a lossless typed superset of every supported format and make each writer return a typed error for values outside its schema.
- O2. Model only the intersection shared by all formats and reject input that cannot enter that common subset. This makes the boundary uniform but excludes valid format-specific data.
- O3. Permit explicit lossy conversion modes whose reports enumerate every dropped value. This can complement either typed representation, but it cannot be the implicit default.

**Recommendation:** O1 as the default contract, with O3 available only through an explicit caller-selected policy. Valid data must never disappear silently.

### A3. Malformed typed fields

- O1. Parse strings into typed values at the input boundary and fail the complete conversion with field and node context.
- O2. Preserve both raw and parsed representations and defer the error to writers that require the typed value. This retains malformed source text but allows invalid state deeper into the application.

**Recommendation:** O1. Inner TreeIR values should be valid by construction; raw source retention belongs in a separate diagnostic structure when required.

Explicit approval is required before recording any architectural decision. Until then, no implementation ticket is ready.

## Locations

- `type TreeIrGraph` [packages/treetime-io/src/graph.rs#L23](../../packages/treetime-io/src/graph.rs#L23)
- `pub fn write_tree_outputs()` [packages/treetime-io/src/graph.rs#L61](../../packages/treetime-io/src/graph.rs#L61)
- `mod tree_ir` adapters [packages/treetime-io/src/tree_ir/mod.rs](../../packages/treetime-io/src/tree_ir/mod.rs)
- `pub fn build_ir_with_mutations()` and sibling projections [packages/treetime/src/commands/shared/ir_projection.rs#L29](../../packages/treetime/src/commands/shared/ir_projection.rs#L29)

## Related KB items

- [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md)
- [kb/proposals/output-format-selection.md](../proposals/output-format-selection.md)
- [kb/proposals/tree-format-model-embedding.md](../proposals/tree-format-model-embedding.md)
- [N-io-auspice-number-formatting-missing-augur-golden-master.md](N-io-auspice-number-formatting-missing-augur-golden-master.md)
- [M-core-mutation-representation-and-format-projection-inconsistent.md](M-core-mutation-representation-and-format-projection-inconsistent.md)
- [M-io-usher-mat-mutation-loss-is-implicit.md](M-io-usher-mat-mutation-loss-is-implicit.md)
- [N-io-phyloxml-mutation-property-contract-undecided.md](N-io-phyloxml-mutation-property-contract-undecided.md)
- [M-timetree-tree-output-inference-metadata-incomplete.md](M-timetree-tree-output-inference-metadata-incomplete.md)
- [M-io-auspice-trait-confidence-silently-coerced.md](M-io-auspice-trait-confidence-silently-coerced.md)
- [M-io-phyloxml-booleans-silently-coerced.md](M-io-phyloxml-booleans-silently-coerced.md)
