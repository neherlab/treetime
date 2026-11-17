# The annotation of ancestral translations is not checked against the alignment

GFF3 coordinates are reference coordinates, and nucleotide mutations are positions in alignment columns. The two agree only when the alignment columns are reference coordinates, with insertions relative to the reference removed.

The partition schemes check this: with `--partition-by codon-positions` or `genes`, each annotation file must place all its CDSs on one sequence, carry a `##sequence-region` for that sequence from 1 to the length of its alignment, and end every CDS segment at or before the last site, or the run stops with an argument error on `annotation` and `alignment` (`fn check_annotation()` in [packages/app-commands/src/commands/shared/partition_by.rs#L458](../../packages/app-commands/src/commands/shared/partition_by.rs#L458)).

The translations of ancestral (`--translations` with `--annotation`) do not run this check. They read the annotation through `fn read_gff3_annotations()` ([packages/app-commands/src/commands/ancestral/aa_node_data.rs#L194](../../packages/app-commands/src/commands/ancestral/aa_node_data.rs#L194)) and check only that each peptide has one third of the CDS length (`fn validate_cds_length()`, [packages/treetime/src/ancestral/aa.rs#L112](../../packages/treetime/src/ancestral/aa.rs#L112)). With an alignment that has insertions relative to the reference, the annotation and the mutations of the node data use different coordinates without a message, and a consumer that maps nucleotide mutations onto the annotation reads them in the wrong codon. A decoder that cuts codons from alignment columns ([M-ancestral-aa-outputs-disagree-with-translated-nucleotide-outputs.md](M-ancestral-aa-outputs-disagree-with-translated-nucleotide-outputs.md)) would shift every codon after the first insertion.

The bundled annotated datasets fit the contract: every bundled annotation carries a `##sequence-region` from 1 to the alignment length.

## Fix direction

Run the check of the schemes on the translation path too, so that both uses of `--annotation` apply one rule.

> [!IMPORTANT]
> **Decision required.** The scheme path requires the `##sequence-region` line. Requiring it for translations stops runs whose annotation has no landmark length and that run today; the options are to require it on both paths or to skip the length check on the translation path when the line is absent.

## Validation

- An alignment with one inserted column and the matching annotation stops a translation run with an error that names both lengths
- The annotated datasets of `data/` pass
