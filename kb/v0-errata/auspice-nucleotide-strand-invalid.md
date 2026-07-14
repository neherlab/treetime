# Auspice nucleotide annotation uses an invalid strand value

## v0 location

`create_auspice_json()` writes `meta.genome_annotations.nuc.strand` as `"+:"` [packages/legacy/treetime/treetime/CLI_io.py#L281-L287](../../packages/legacy/treetime/treetime/CLI_io.py#L281-L287).

## Erratum

The value contains an extra colon. It is outside the Augur annotation schema and does not express any additional coordinate information.

## Evidence

- Augur's annotation schema allows only `"+"` for the optional nucleotide strand and requires the nucleotide interval to begin at one [[src](https://github.com/nextstrain/augur/blob/d8e38736037ba9474a809f9a5a63bc2b279d2407/augur/data/schema-annotations.json#L7-L24)].
- Auspice treats the nucleotide annotation as a positive-strand genome unless its strand is exactly `"-"`; the strand does not encode frame or delimiter information [[src](https://github.com/nextstrain/auspice/blob/0bb35258e105932a1f2d3839049f13292ff2d14b/src/util/entropyCreateStateFromJsons.ts#L108-L125)].
- V0's neighboring `type: "source"`, `start: 1`, and `end: full_length` fields already describe the annotation independently of the malformed strand string [packages/legacy/treetime/treetime/CLI_io.py#L281-L287](../../packages/legacy/treetime/treetime/CLI_io.py#L281-L287).

## v0 impact

- The generated document fails strict validation against the Augur annotation schema.
- Current Auspice accepts the nucleotide track because it rejects only exact negative strand, but other schema-aware consumers may reject the document.

## v1 status

V1 does not yet emit nucleotide genome annotations. [kb/issues/M-timetree-tree-output-inference-metadata-incomplete.md](../issues/M-timetree-tree-output-inference-metadata-incomplete.md) and [kb/tickets/timetree-complete-tree-output-inference-metadata.md](../tickets/timetree-complete-tree-output-inference-metadata.md) require schema-valid `"+"` output.
