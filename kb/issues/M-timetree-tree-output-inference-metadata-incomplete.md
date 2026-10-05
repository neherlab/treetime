# TreeTime tree-output inference metadata is incomplete

> [!WARNING]
> **Needs review.** The `node_attrs.confidence` that v1 writes is the input tree's branch support, not the v0 mutation-derived value. `fn with_branch_support()` [`packages/app-output/src/tree_output.rs#L201-L221`](../../packages/app-output/src/tree_output.rs#L201-L221) writes `branch_support`, which `fn timetree_node_outputs()` takes from the input confidences ([`packages/app-commands/src/commands/timetree/run.rs#L50`](../../packages/app-commands/src/commands/timetree/run.rs#L50), [`#L407`](../../packages/app-commands/src/commands/timetree/run.rs#L407)). No code in `packages/app-output/src` computes $1 - e^{-m}$, and the timetree colorings [`packages/app-output/src/timetree_tree_output.rs#L69-L70`](../../packages/app-output/src/timetree_tree_output.rs#L69-L70) contain no `confidence` key. Confirm whether input support should be kept, replaced, or emitted beside the v0 value.

TreeTime writes Auspice directly from the graph in [`packages/app-output/src/tree_output.rs`](../../packages/app-output/src/tree_output.rs), but the Auspice metadata payload is not yet a complete output contract.

Present today: every command with a root sequence writes `meta.genome_annotations.nuc` (start 1, end at the alignment length, strand `"+"`, type `"source"`) and declares the `entropy` panel ([`packages/app-output/src/tree_output.rs#L130-L134`](../../packages/app-output/src/tree_output.rs#L130-L134)). When the input tree carries branch support, it is written as `node_attrs.confidence`.

## Problem

- Insertion and deletion handling across formats is governed by [M-core-mutation-representation-and-format-projection-inconsistent.md](M-core-mutation-representation-and-format-projection-inconsistent.md).
- The mutation-derived branch support below is not computed; nodes of an input tree without support values get no `node_attrs.confidence`.
- `meta.colorings` omits the continuous `confidence` coloring ("Branch Support") that v0 declares.

## Branch support formula

For $m$ canonical nucleotide substitutions (`A`, `C`, `G`, `T`) on a branch, v0 reports confidence

$$c = 1 - e^{-m}$$

where $c$ is branch support confidence, rounded to three decimals. Leaves receive $c=1$. v0 writes it only when run with sequence data, and declares the coloring `{'title': 'Branch Support', 'type': 'continuous', 'key': 'confidence'}` ([`packages/legacy/treetime/treetime/CLI_io.py#L300`](../../packages/legacy/treetime/treetime/CLI_io.py#L300), [`packages/legacy/treetime/treetime/CLI_io.py#L318-L327`](../../packages/legacy/treetime/treetime/CLI_io.py#L318-L327)).

## Potential solutions

- O1. Extend the shared projection with all mutation and inference metadata needed by Auspice.
- O2. Build format-specific TreeTime projections. This can preserve format-specific data but duplicates the shared boundary.

## Recommendation

Extend the TreeTime projection once, supplying typed Auspice mutations, confidence, alignment length, annotations, and numerical dates with confidence intervals. Add an Auspice fixture for the complete visualization payload.

## Fix

1. Consume the canonical typed nucleotide mutation projection for Auspice output.
2. Compute branch support as $c=1-e^{-m}$, where $m$ is the number of canonical nucleotide substitutions; set leaf confidence to $1$.
3. Emit `node_attrs.confidence` and `meta.colorings.confidence`.
4. Pass alignment length and annotation data through the command projection rather than recovering them in the writer.
5. Keep confidence, annotations, panels, and mutation output mutually consistent when an optional source is absent.

## Validation

- Whole-document golden master against pinned Augur/Auspice behavior
- Internal-node and leaf confidence cases
- Auspice substitution, insertion, and deletion projection cases for verified mappings
- Alignment with and without genome annotations

## Related issues

- [M-core-mutation-representation-and-format-projection-inconsistent.md](M-core-mutation-representation-and-format-projection-inconsistent.md)
- [N-ancestral-auspice-json-not-produced.md](N-ancestral-auspice-json-not-produced.md)
- [kb/reports/augur-node-data-json.md](../reports/augur-node-data-json.md)
