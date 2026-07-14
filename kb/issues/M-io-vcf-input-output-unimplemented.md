# VCF input and output are unimplemented

TreeTime v1 exposes `--vcf-reference` on sequence-consuming commands, but it has no VCF reader or writer. VCF input, compressed VCF input, and VCF output remain unchecked in [kb/features/io.md](../features/io.md).

## User-facing impact

Variant-only workflows cannot pass a VCF plus reference sequence directly to v1 or receive reconstructed variants as VCF. Users must materialize a full alignment before analysis and convert output back to variants with another tool. This loses the memory advantage of a sparse variant representation for large genomes and sample collections.

Public reports describe two related workflows:

- A SNP-only alignment user asks whether timetree analysis requires VCF input or manual rate rescaling [[issue](https://github.com/neherlab/treetime/issues/316)].
- A Python TreeTime user reports that a multi-record reference FASTA is rejected for a multi-chromosome VCF workflow [[issue](https://github.com/neherlab/treetime/issues/247)]. A maintainer confirms that the Python implementation does not handle multi-chromosome VCF files [[comment](https://github.com/neherlab/treetime/issues/247#issuecomment-1625603477)]. This is related format evidence, not proof of a shared implementation defect in Rust.

## Scope axes requiring decisions

### Parser and compression

Select a maintained VCF parser and define support for plain and compressed VCF streams, including headers, genotype fields, symbolic alleles, missing calls, and non-UTF-8 paths.

### Internal representation

Define whether VCF records become a sparse partition directly or are expanded into aligned sequences. The representation must preserve reference coordinates, missing and ambiguous calls, insertions, deletions, and per-sample variants without unnecessary dense materialization.

### Output contract

Define which reconstructed states are emitted, how root and ancestral samples are represented, how unchanged and missing positions are handled, and whether output preserves input headers and contig metadata.

### Multi-contig semantics

Decide whether contigs are independent partitions, one coordinate space with contig-qualified positions, or separate command invocations. The choice affects reference loading, sample identity reconciliation, output ordering, and model assignment.

## Ticket readiness

No implementation ticket is ready. The parser, compression, internal representation, output contract, and multi-contig semantics determine incompatible APIs and data models. A ticket can be created after these axes are decided and the accepted design specifies end-to-end behavior and validation oracles.

## Related

- [kb/features/io.md](../features/io.md) - VCF input, compressed VCF input, and VCF output inventory
- [M-clock-dead-cli-arguments.md](M-clock-dead-cli-arguments.md) - `--vcf-reference` is parsed but unused by `clock`
- [N-timetree-dead-cli-flags.md](N-timetree-dead-cli-flags.md) - `--vcf-reference` is parsed but unused by `timetree`
- [N-io-large-dataset-memory-constraint.md](N-io-large-dataset-memory-constraint.md) - dense alignment materialization increases peak memory
- [N-io-multi-segment-genome-input.md](N-io-multi-segment-genome-input.md) - partition semantics for segmented genomes
