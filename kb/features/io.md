# I/O Formats

## Input Formats

- [x] FASTA (plain text)
- [x] Compressed FASTA (gzip, xz, bzip2, zstd)
- [x] Several FASTA files, numbered in the order of the `--alignment` list, one partition per file unless `--partition` entries split the files into sets of sites ([kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md)); homoplasy, `--translations` and `--annotation` accept one file
- [/] Stdin input (`-` for at most one input of a run; `--config -` together with an input `-` reads standard input twice: [kb/issues/M-io-stdin-config-and-input-read-twice.md](../issues/M-io-stdin-config-and-input-read-twice.md))
- [x] Newick tree
- [x] Nexus tree
- [ ] PHYLIP alignment ([kb/issues/M-io-phylip-alignment-input-unsupported.md](../issues/M-io-phylip-alignment-input-unsupported.md))
- [x] CSV/TSV dates
- [x] CSV/TSV discrete states (mugration)
- [ ] VCF input (variant call format)
- [ ] Compressed VCF (.vcf.gz)
- [ ] Custom GTR model file

## Output Formats

- [x] FASTA (ancestral sequences)
- [x] Newick (plain, or with BEAST or NHX node comments through `--newick-annotations`)
- [x] Nexus (always with BEAST node comments: mutations, date, trait)
- [x] JSON (model document, clock model, graph)
- [x] Per-file sequence outputs of several alignment files (FASTA, Auspice JSON, UShER MAT, Nexus, Newick with `--newick-annotations beast` or `nhx`, ancestral augur node data): `<command>.<n>.<ext>` under `--output-all`, or a `{alignment}` placeholder in the per-file flag; a per-file flag without the placeholder stops with an error
- [x] CSV (clock regression, confidence intervals)
- [x] SVG/PNG/PDF charts (clock regression)
- [x] Graphviz DOT
- [x] Output topology ordering
- [ ] VCF output (v0 writes .vcf for VCF inputs: [kb/issues/M-io-vcf-input-output-unimplemented.md](../issues/M-io-vcf-input-output-unimplemented.md))
- [/] Auspice v2 JSON (schema-validated for all tree-writing commands; no Augur golden master checks the document: [kb/issues/N-io-auspice-output-missing-augur-golden-master.md](../issues/N-io-auspice-output-missing-augur-golden-master.md))
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

- [x] `--config` per-command configuration object (JSON or YAML, optionally compressed), in the shape the command serializes to. Precedence: explicit CLI flag, then config file, then default. Scope is a single command's own args. The `partition` key of ancestral, optimize and timetree holds the partition entries (sites and model of each partition, [kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md)); a config file for several commands or for partitions of other data types remains an open design in [kb/proposals/config-file-multi-partition.md](../proposals/config-file-multi-partition.md).
