# TreeTime tree-output inference metadata is incomplete

TreeTime writes Auspice and PhyloXML directly from the graph in [`packages/treetime/src/commands/shared/tree_output.rs`](../../packages/treetime/src/commands/shared/tree_output.rs), but the metadata payload is not yet a complete output contract.

> [!NOTE]
> The tree-output refactor added several of the fields this issue tracked: `treetime:gamma`, `treetime:input_branch_support`/`confidence`, `date_is_inferred`, and `genome_annotations` are now emitted. The **remaining** items below are not confirmed present; re-audit each against `tree_output.rs`.

## Problem

- Insertion and deletion handling across formats is governed by [M-core-mutation-representation-and-format-projection-inconsistent.md](M-core-mutation-representation-and-format-projection-inconsistent.md).
- `meta.colorings` may still omit the continuous `confidence` coloring that corresponds to branch support.
- `meta.genome_annotations.nuc` completeness (alignment coordinates, type, strand) is unconfirmed.
- `date_raw_value` is not emitted; confirm whether `date_is_inferred` and relaxed-clock `gamma` are populated from command state in every mode rather than left at defaults.
- PhyloXML inferred-date and `gamma` properties for relaxed-clock output are unconfirmed.

The Auspice fields share one command projection and must be validated together. The already-defined non-mutation PhyloXML date and gamma properties use the same command state. PhyloXML mutation vocabulary and indel representation remain excluded pending [N-io-phyloxml-mutation-property-contract-undecided.md](N-io-phyloxml-mutation-property-contract-undecided.md).

## Branch support formula

For $m$ canonical nucleotide substitutions on a branch, v0 reports confidence

$$c = 1 - e^{-m}$$

where $c$ is branch support confidence. Leaves receive $c=1$.

## Potential solutions

- O1. Extend the shared TreeIR projection with all mutation and inference metadata needed by Auspice and PhyloXML.
- O2. Build format-specific TreeTime projections. This can preserve format-specific data but duplicates the shared boundary.

## Recommendation

Extend the TreeTime projection once, supplying typed Auspice mutations, confidence, alignment length, annotations, raw/inferred dates, and relaxed-clock rate multipliers. Add an Auspice fixture for the complete visualization payload and a PhyloXML fixture limited to the already-defined non-mutation inference properties.

## Related issues

- [M-core-mutation-representation-and-format-projection-inconsistent.md](M-core-mutation-representation-and-format-projection-inconsistent.md)
- [N-ancestral-auspice-json-not-produced.md](N-ancestral-auspice-json-not-produced.md)
- [N-io-phyloxml-mutation-property-contract-undecided.md](N-io-phyloxml-mutation-property-contract-undecided.md)
- [kb/reports/augur-node-data-json.md](../reports/augur-node-data-json.md)
