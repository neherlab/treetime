# TreeTime tree-output inference metadata is incomplete

TreeTime writes Auspice directly from the graph in [`packages/app-output/src/auspice.rs`](../../packages/app-output/src/auspice.rs), but the Auspice metadata payload is not yet a complete output contract.

Present today: every command with a root sequence writes `meta.genome_annotations.nuc` (start 1, end at the alignment length, strand `"+"`, type `"source"`) and declares the `entropy` panel ([`packages/app-output/src/auspice.rs#L140-L144`](../../packages/app-output/src/auspice.rs#L140-L144)).

Branch support, both the input values and the mutation-based value that v0 writes as `confidence`, is tracked in [M-io-branch-support-dropped-from-outputs.md](M-io-branch-support-dropped-from-outputs.md).

## Problem

- Insertion and deletion handling across formats is governed by [M-core-mutation-representation-and-format-projection-inconsistent.md](M-core-mutation-representation-and-format-projection-inconsistent.md).

## Potential solutions

- O1. Extend the shared projection with all mutation and inference metadata needed by Auspice.
- O2. Build format-specific TreeTime projections. This can preserve format-specific data but duplicates the shared boundary.

## Recommendation

Extend the TreeTime projection once, supplying typed Auspice mutations, alignment length, annotations, and numerical dates with confidence intervals. Add an Auspice fixture for the complete visualization payload.

## Fix

1. Consume the canonical typed nucleotide mutation projection for Auspice output.
2. Pass alignment length and annotation data through the command projection rather than recovering them in the writer.
3. Keep annotations, panels, and mutation output mutually consistent when an optional source is absent.

## Validation

- Whole-document golden master against pinned Augur/Auspice behavior
- Auspice substitution, insertion, and deletion projection cases for verified mappings
- Alignment with and without genome annotations

## Related issues

- [M-core-mutation-representation-and-format-projection-inconsistent.md](M-core-mutation-representation-and-format-projection-inconsistent.md)
- [N-ancestral-auspice-json-not-produced.md](N-ancestral-auspice-json-not-produced.md)
- [kb/reports/augur-node-data-json.md](../reports/augur-node-data-json.md)
