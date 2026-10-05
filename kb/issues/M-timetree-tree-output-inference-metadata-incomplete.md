# TreeTime tree-output inference metadata is incomplete

> [!WARNING]
> **Needs review.** v1 writes no mutation-derived `node_attrs.confidence`. `fn node_attrs_other()` [`packages/app-output/src/auspice.rs#L337-L357`](../../packages/app-output/src/auspice.rs#L337-L357) writes the input tree's branch support as `node_attrs.confidence` when the struct of facts carries `branch_support`, and timetree carries none ([`packages/app-commands/src/commands/timetree/run.rs#L349`](../../packages/app-commands/src/commands/timetree/run.rs#L349)), because support values sit on the wrong splits after a reroot ([M-io-input-branch-support-not-moved-on-reroot.md](M-io-input-branch-support-not-moved-on-reroot.md)). No code in `packages/app-output/src` computes $1 - e^{-m}$, and the colorings [`packages/app-output/src/auspice.rs#L116-L128`](../../packages/app-output/src/auspice.rs#L116-L128) contain no `confidence` key. Confirm whether input support should be kept, replaced, or emitted beside the v0 value.

TreeTime writes Auspice directly from the graph in [`packages/app-output/src/auspice.rs`](../../packages/app-output/src/auspice.rs), but the Auspice metadata payload is not yet a complete output contract.

Present today: every command with a root sequence writes `meta.genome_annotations.nuc` (start 1, end at the alignment length, strand `"+"`, type `"source"`) and declares the `entropy` panel ([`packages/app-output/src/auspice.rs#L140-L144`](../../packages/app-output/src/auspice.rs#L140-L144)). Timetree writes no input branch support.

## Problem

- Insertion and deletion handling across formats is governed by [M-core-mutation-representation-and-format-projection-inconsistent.md](M-core-mutation-representation-and-format-projection-inconsistent.md).
- The mutation-derived branch support below is not computed, so timetree nodes get no `node_attrs.confidence`.
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
