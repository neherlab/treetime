# TreeTime tree-output inference metadata is incomplete

TreeTime writes Auspice directly from the graph in [`packages/treetime/src/commands/shared/tree_output.rs`](../../packages/treetime/src/commands/shared/tree_output.rs), but the Auspice metadata payload is not yet a complete output contract.

> [!NOTE]
> The tree-output refactor added several of the fields this issue tracked: branch-support `confidence` and `genome_annotations` are now emitted. The **remaining** items below are not confirmed present; re-audit each against `tree_output.rs`.

## Problem

- Insertion and deletion handling across formats is governed by [M-core-mutation-representation-and-format-projection-inconsistent.md](M-core-mutation-representation-and-format-projection-inconsistent.md).
- `meta.colorings` may still omit the continuous `confidence` coloring that corresponds to branch support.
- `meta.genome_annotations.nuc` completeness (alignment coordinates, type, strand) is unconfirmed.

## Branch support formula

For $m$ canonical nucleotide substitutions on a branch, v0 reports confidence

$$c = 1 - e^{-m}$$

where $c$ is branch support confidence. Leaves receive $c=1$.

## Potential solutions

- O1. Extend the shared projection with all mutation and inference metadata needed by Auspice.
- O2. Build format-specific TreeTime projections. This can preserve format-specific data but duplicates the shared boundary.

## Recommendation

Extend the TreeTime projection once, supplying typed Auspice mutations, confidence, alignment length, annotations, and numerical dates with confidence intervals. Add an Auspice fixture for the complete visualization payload.

## Related issues

- [M-core-mutation-representation-and-format-projection-inconsistent.md](M-core-mutation-representation-and-format-projection-inconsistent.md)
- [N-ancestral-auspice-json-not-produced.md](N-ancestral-auspice-json-not-produced.md)
- [kb/reports/augur-node-data-json.md](../reports/augur-node-data-json.md)
