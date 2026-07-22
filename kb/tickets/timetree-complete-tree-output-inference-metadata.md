# Complete TreeTime tree-output inference metadata

Produce internally consistent TreeTime Auspice metadata and populate the already-defined non-mutation PhyloXML inference properties.

> [!NOTE]
> This ticket predates the tree-output refactor, which removed the `tree_ir` layer and now writes formats directly from the graph in `packages/treetime/src/commands/shared/tree_output.rs` (readers in `packages/treetime-io/src/auspice.rs` and `packages/treetime-io/src/phyloxml.rs`). Its parent issue's current status is unconfirmed; re-derive the steps and code locations against the current code before executing.

## Required changes

1. Consume the canonical typed nucleotide mutation projection for Auspice output.
2. Compute branch support as $c=1-e^{-m}$, where $m$ is the number of canonical nucleotide substitutions; set leaf confidence to $1$.
3. Emit `node_attrs.confidence`, `meta.colorings.confidence`, `meta.genome_annotations.nuc`, and the corresponding panel metadata. The nucleotide annotation must use schema-valid strand `"+"`, not v0's malformed `"+:"` value.
4. Pass alignment length and annotation data through the command projection rather than recovering them in the writer.
5. Keep confidence, annotations, panels, and mutation output mutually consistent when an optional source is absent.
6. Populate `date_is_inferred`, `date_raw_value`, and relaxed-clock `gamma` from existing TreeTime command state.
7. Populate the existing PhyloXML `date_is_inferred`, `date_raw_value`, and relaxed-clock `gamma` properties. Do not change PhyloXML mutation or indel properties until their vocabulary is approved.

## Validation

- Whole-document golden master against pinned Augur/Auspice behavior.
- Internal-node and leaf confidence cases.
- Auspice substitution, insertion, and deletion projection cases for verified mappings.
- Alignment with and without genome annotations.
- Raw date, inferred date, and relaxed-clock gamma cases in PhyloXML.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/M-timetree-tree-output-inference-metadata-incomplete.md](../issues/M-timetree-tree-output-inference-metadata-incomplete.md)
- Related: [kb/issues/N-io-phyloxml-mutation-property-contract-undecided.md](../issues/N-io-phyloxml-mutation-property-contract-undecided.md)
