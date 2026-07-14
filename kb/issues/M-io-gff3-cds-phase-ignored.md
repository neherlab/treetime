# GFF3 CDS phase is ignored and not validated

The GFF3 reader requires nine tab-separated columns but never parses the CDS phase column. It reads coordinates, strand, and attributes, then constructs `RawCdsRow` without phase data. [`packages/treetime-io/src/gff.rs#L75-L131`](../../packages/treetime-io/src/gff.rs#L75-L131)

`GffCdsSegment` retains only start and end coordinates, so phase cannot be preserved when segmented CDS annotations are emitted. [`packages/treetime-io/src/gff.rs#L37-L42`](../../packages/treetime-io/src/gff.rs#L37-L42) [`packages/treetime-io/src/gff.rs#L162-L168`](../../packages/treetime-io/src/gff.rs#L162-L168)

Malformed phase values are accepted as valid GFF3, and valid nonzero phases are silently discarded. The GFF3 specification requires every CDS phase to be `0`, `1`, or `2` and defines it relative to the CDS feature's 5′ end [[spec](https://github.com/The-Sequence-Ontology/Specifications/blob/fe73505276dd324bf6a55773f3413fe2bed47af4/gff3.md#L59-L62)]. For split CDS features, losing phase makes the annotation insufficient to reconstruct the original CDS coordinate semantics.

## Decision required

Validation can independently reject values outside the GFF3 phase domain. Preserving valid phase requires confirming how the public Augur node-data annotation schema represents per-segment phase. No implementation ticket is ready until that output contract is established.
