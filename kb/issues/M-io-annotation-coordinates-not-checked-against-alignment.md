# Annotation coordinates are not checked against the alignment

GFF3 coordinates are reference coordinates, and nucleotide mutations are positions in alignment columns. The two agree only when the alignment columns are reference coordinates, with insertions relative to the reference removed. TreeTime does not check this:

- The GFF3 reader checks only that each row has $1 \le \text{start} \le \text{end}$ and that the CDS length is a multiple of 3 ([packages/treetime-io/src/gff.rs#L40-L136](../../packages/treetime-io/src/gff.rs#L40-L136))
- The amino-acid route checks only that each peptide has one third of the CDS length (`fn validate_cds_length()` [packages/treetime/src/ancestral/aa.rs#L112-L127](../../packages/treetime/src/ancestral/aa.rs#L112-L127))
- The node data carries the annotation coordinates next to nucleotide mutations in alignment columns

With an alignment that has insertions relative to the reference, the annotation and the mutations of the node data use different coordinates without a message, and a consumer that maps nucleotide mutations onto the annotation reads them in the wrong codon. A decoder that cuts codons from alignment columns ([M-ancestral-aa-outputs-disagree-with-translated-nucleotide-outputs.md](M-ancestral-aa-outputs-disagree-with-translated-nucleotide-outputs.md)) would shift every codon after the first insertion.

The bundled annotated datasets fit the contract: on `rsv/a/2000` and `mpox/clade-ii/2000` the alignment length equals the reference length (15,225 and 197,209), and every CDS is one row with phase 0.

## Fix direction

Compare the alignment length with the length of the landmark sequence (`##sequence-region` or the `region` row of the annotation) and stop with an error when they differ. Document that alignment columns must be reference coordinates. Axis 6 of [kb/proposals/codon-aware-ancestral-reconstruction.md](../proposals/codon-aware-ancestral-reconstruction.md) describes the whole reader.

> [!IMPORTANT]
> **Decision required.** A length mismatch can stop the run with an error or continue with a warning in the run record. An error blocks alignments that carry extra columns at the ends without insertions in coding regions; a warning lets shifted codons through. An annotation without a landmark length allows no check at all; the options are to require the landmark length or to skip the check with a warning.

## Validation

- An alignment with one inserted column and the matching annotation stops with an error that names both lengths
- The annotated datasets of `data/` pass
