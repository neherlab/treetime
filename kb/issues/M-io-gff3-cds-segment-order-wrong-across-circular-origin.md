# GFF3 CDS segments are out of order for a CDS across the origin of a circular genome

The GFF3 reader joins the `CDS` rows that share an `ID` into one CDS, sorts its segments by coordinate, and reverses them for the minus strand (`fn assemble_feature()`, [packages/treetime-io/src/gff.rs#L280-L301](../../packages/treetime-io/src/gff.rs#L280-L301)). This gives the biological 5' to 3' order for every CDS that does not wrap.

A CDS that crosses the origin of a circular genome (a plasmid, a mitochondrion, a circular virus) has a segment at the end of the sequence and a segment at its start. GFF3 writes such a CDS either as two rows that end at the last site and start at site 1, or as one row whose end exceeds the landmark length ([GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/fe73505276dd324bf6a55773f3413fe2bed47af4/gff3.md?plain=1), "Circular genomes"). For the two-row form, the sort puts the segment at site 1 first on the forward strand, so the codons of the CDS start in the wrong segment. Nextclade keeps the segments in file order and wraps circular genomes ([cds.rs#L46-L90](https://github.com/nextstrain/nextclade/blob/6190424cf8014d999fccf59d049c2c78a8bc28bb/packages/nextclade/src/gene/cds.rs#L46-L90)).

The partition schemes (`--partition-by codon-positions`) take the codon position of each site from the phase of its own segment, so their codon positions stay correct for any segment order. The translations of ancestral (`--translations`) and the segments of the augur node data annotation read the segments in the stored order.

> [!IMPORTANT]
> **Investigation required.** No bundled dataset has a circular genome, so the two-row form is not reproduced in TreeTime. Check a circular annotation (for example a plasmid from RefSeq) against the translations of ancestral, and decide whether the reader keeps file order, detects the wrap from the landmark length, or rejects a CDS across the origin.

## Fix direction

Detect the wrap from the `Is_circular=true` attribute of the landmark or from a segment that ends at the landmark end next to a segment that starts at site 1, and order the segments 5' to 3' across the origin. A one-row CDS with an end past the landmark length needs the same treatment, or a clear error.

## Validation

- A forward-strand CDS of two rows across the origin gives its segments in 5' to 3' order, and its translation reads the codon that spans the origin
- The annotations of `data/` give the same CDSs as before

## Related

- [M-ancestral-node-data-annotation-drops-cds-phase.md](M-ancestral-node-data-annotation-drops-cds-phase.md): another property of CDS segments in the same outputs

## Proposed implementation plan

The investigation above decides between three reader behaviors. Steps 1 and 3 hold for each of them; step 2 has one variant per option.

1. Data: the reader already keeps `##sequence-region` directives (`GffSequenceRegion`, `fn parse_sequence_region()` in `packages/treetime-io/src/gff.rs`). Read `Is_circular=true` from the landmark `region` row as well, so that `fn assemble_feature()` knows the landmark length and circularity of `seqid`
2. Ordering in `fn assemble_feature()` (`packages/treetime-io/src/gff.rs`), after the coordinate sort:
   - Wrap detection: on a circular landmark, when a segment ends at the landmark end and another starts at site 1, rotate the sorted segments so that the segment after the largest gap comes first, then reverse for the minus strand. Split a one-row CDS with `end > landmark length` into two segments the same way
   - File order: keep rows in `line_no` order on the forward strand and reversed order on the minus strand (Nextclade behavior)
   - Reject: return a `make_error!` that names the CDS and the origin when a CDS on a circular landmark touches both ends, or ends past the landmark length
3. Tests in `packages/treetime-io/src/__tests__/test_gff.rs`, inline GFF3 literals:
   - Forward-strand CDS of two rows `91..100` and `1..5` on a 100-site circular landmark: expected segments `[91..100, 1..5]` (oracle: GFF3 specification, "Circular genomes")
   - Same CDS on the minus strand: segments `[1..5, 91..100]`
   - One-row CDS `91..105` on a 100-site landmark: two segments or the error, per the chosen option
   - Linear landmark with the same rows keeps the coordinate order (no wrap without `Is_circular`)
4. Validate: `just lint-rs`, `just test-rs -E 'package(treetime-io)'`, and `just smoke` on a dataset with an annotation to confirm that the CDSs of `data/` are unchanged
