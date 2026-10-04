# I/O Formats

## Input Formats

- [x] FASTA (plain text)
- [x] Compressed FASTA (gzip, xz, bzip2, zstd)
- [x] Multiple FASTA files (concatenated)
- [x] Stdin input
- [x] Newick tree
- [ ] Nexus tree ([kb/issues/M-io-nexus-and-phylip-input-unsupported.md](../issues/M-io-nexus-and-phylip-input-unsupported.md))
- [ ] PHYLIP alignment ([kb/issues/M-io-nexus-and-phylip-input-unsupported.md](../issues/M-io-nexus-and-phylip-input-unsupported.md))
- [x] CSV/TSV dates
- [x] CSV/TSV discrete states (mugration)
- [ ] VCF input (variant call format)
- [ ] Compressed VCF (.vcf.gz)
- [ ] Custom GTR model file

## Output Formats

- [x] FASTA (ancestral sequences)
- [x] Newick (annotated tree)
- [x] Nexus (annotated tree)
- [x] JSON (GTR model, clock model, graph)
- [x] CSV (clock regression, confidence intervals)
- [x] SVG/PNG charts (clock regression)
- [x] Graphviz DOT
- [x] Output topology ordering
- [ ] VCF output (v0 writes .vcf for VCF inputs: [kb/issues/M-io-vcf-input-output-unimplemented.md](../issues/M-io-vcf-input-output-unimplemented.md))
- [/] Auspice v2 JSON (schema-validated for all tree-writing commands; entropy perturbs the Shannon definition and inference metadata is incomplete: [kb/issues/M-io-auspice-entropy-perturbs-shannon-definition.md](../issues/M-io-auspice-entropy-perturbs-shannon-definition.md), [kb/issues/M-timetree-tree-output-inference-metadata-incomplete.md](../issues/M-timetree-tree-output-inference-metadata-incomplete.md))
- [ ] Skyline TSV/plot
- [ ] Substitution rates TSV
- [ ] Outliers TSV
- [ ] Branch mutations TXT
- [ ] Sequence evolution model TXT

## v1-Only Formats

- [/] PhyloXML (format type model, reader, and writer in `util-phyloxml`; no command reads or writes PhyloXML: [kb/issues/N-io-phyloxml-crate-unused.md](../issues/N-io-phyloxml-crate-unused.md))
- [x] UShER MAT output (validates global reference nucleotides and writes indels as missing data: [kb/decisions/io-usher-mat-gaps-as-missing-data.md](../decisions/io-usher-mat-gaps-as-missing-data.md); holds the nucleotide mutations only, because MAT states are A, C, G and T [[src](https://github.com/yatisht/usher/blob/ac9c982d937c3bc43e2c16a0b73d30cf2b937118/parsimony.proto#L4-L11)], and UShER tools compute amino-acid changes from them with a GTF annotation)
- [ ] UShER MAT input ([kb/issues/N-io-usher-mat-input-not-wired.md](../issues/N-io-usher-mat-input-not-wired.md))
- [ ] Taxonium JSONL output ([kb/issues/N-io-taxonium-jsonl-output-unsupported.md](../issues/N-io-taxonium-jsonl-output-unsupported.md))
- [ ] MAPLE alignment input ([kb/issues/N-io-maple-alignment-input-unsupported.md](../issues/N-io-maple-alignment-input-unsupported.md))
- [ ] Nextclade NDJSON alignment input ([kb/issues/N-io-nextclade-ndjson-input-unsupported.md](../issues/N-io-nextclade-ndjson-input-unsupported.md))
- [ ] tskit tree sequence input ([kb/issues/N-io-tskit-tree-sequence-input-unsupported.md](../issues/N-io-tskit-tree-sequence-input-unsupported.md))
- [x] YAML serialization
- [x] Compressed FASTA output
- [x] Streaming readers/writers with automatic decompression

## Configuration Input

- [x] `--config` per-command configuration object (JSON or YAML, optionally compressed), in the shape the command serializes to. Precedence: explicit CLI flag, then config file, then default. Scope is a single command's own args; the broader multi-partition config system (multiple alignments, per-partition models) remains an open design in [kb/proposals/config-file-multi-partition.md](../proposals/config-file-multi-partition.md).
