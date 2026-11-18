# The augur node data annotation of ancestral drops the CDS phase

The GFF3 reader validates the phase of every `CDS` row as 0, 1 or 2 and stores it on the segment (`struct GffCdsSegment`, [packages/treetime-io/src/gff.rs#L89-L93](../../packages/treetime-io/src/gff.rs#L89-L93)); any other value stops the run with an error that names the value, the CDS and the line. The codon positions of `--partition-by codon-positions` use the phase.

The annotation that ancestral writes into the augur node data (`fn gff_cds_to_annotation()`, [packages/app-commands/src/commands/ancestral/aa_node_data.rs#L240-L277](../../packages/app-commands/src/commands/ancestral/aa_node_data.rs#L240-L277)) writes start, end and strand of each CDS or segment, and no phase. The GFF3 specification defines the phase relative to the 5' end of each row ([spec](https://github.com/The-Sequence-Ontology/Specifications/blob/fe73505276dd324bf6a55773f3413fe2bed47af4/gff3.md#L59-L62)), so a consumer of a split CDS with a nonzero phase cannot rebuild the codons from the node data alone. Every CDS of the bundled annotations has phase 0.

> [!IMPORTANT]
> **Investigation required.** Confirm how the public augur node data schema represents the phase of a CDS segment (`AugurNodeDataJsonAnnotationSegment` in `packages/util-augur-node-data-json/src/annotations.rs` keeps unknown keys in `other`) before choosing the output field.

## Fix direction

Write the phase of each segment into the annotation in the representation that augur and Auspice read, once the schema is confirmed.

## Validation

- A two-segment CDS whose second segment has phase 1 gives node data whose annotation carries phase 1 on that segment
- The annotations of `data/` give the same node data as before

## Related

- [M-io-gff3-cds-segment-order-wrong-across-circular-origin.md](M-io-gff3-cds-segment-order-wrong-across-circular-origin.md): the order of the same segments
