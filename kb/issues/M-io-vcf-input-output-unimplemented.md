# VCF input and output are unimplemented

> [!IMPORTANT]
> **Decision required.** The parser and compression support, internal representation, output contract, and multi-contig semantics are undecided, and they determine incompatible APIs and data models. See "Scope axes requiring decisions".

TreeTime v1 has no VCF reader or writer, and no command has a VCF flag. The 0.x flag `--vcf-reference` stops with a message that it is unavailable (`packages/app-cli/src/cli/v0_args.rs`). VCF input, compressed VCF input, and VCF output remain unchecked in [kb/features/io.md](../features/io.md).

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

## Readiness

The parser, compression, internal representation, output contract, and multi-contig semantics determine incompatible APIs and data models. Start implementation after these axes are decided and the accepted design specifies end-to-end behavior and validation oracles.

## Related

- [kb/features/io.md](../features/io.md) - VCF input, compressed VCF input, and VCF output inventory
- [N-io-large-dataset-memory-constraint.md](N-io-large-dataset-memory-constraint.md) - dense alignment materialization increases peak memory
- [N-io-multi-segment-genome-input.md](N-io-multi-segment-genome-input.md) - partition semantics for segmented genomes

## Proposed implementation plan

The steps stop at the open decision above; each step uses the choice of its axis.

1. Parser: vet the chosen parser crate with the `supply-chain` skill and add it through "Adding a dependency"; add a reader in a new `packages/treetime-io/src/vcf.rs` for plain and gzip streams that takes `PathBuf` input
2. Representation: convert records into the chosen form (sparse partition from reference plus per-sample variants, or aligned sequences), with missing calls, ambiguous calls, and indels kept as decided
3. CLI: add the VCF and reference inputs to the commands that read alignments, and change the `--vcf-reference` entry in `packages/app-cli/src/cli/v0_args.rs` from unavailable to the v1 replacement; regenerate with `just gen`
4. Output: add a VCF writer for reconstructed states per the output contract
5. Tests: reader unit tests on small hand-written VCF text (oracle: the VCF specification); sparse result of VCF plus reference equals the result of the equivalent FASTA alignment (oracle: dense vs sparse comparison); ancestral VCF output on the v0 VCF example data compared with v0 output (oracle: v0, dense)
6. Validate with `just lint-rs`, `just test-rs`, `just generated-check`, and `just r treetime ancestral` on a VCF dataset
