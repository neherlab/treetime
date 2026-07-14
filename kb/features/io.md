# I/O Formats

## Input Formats

- [x] FASTA (plain text)
- [x] Compressed FASTA (gzip, xz, bzip2, zstd)
- [x] Multiple FASTA files (concatenated)
- [x] Stdin input
- [x] Newick tree
- [x] Nexus tree
- [x] Phylip tree
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
- [/] Output topology ordering (direct formats honor the plan; TreeIR-backed formats can project the unordered graph: [kb/issues/M-io-tree-backed-output-order-inconsistent.md](../issues/M-io-tree-backed-output-order-inconsistent.md))
- [ ] VCF output (v0 writes .vcf for VCF inputs: [kb/issues/M-io-vcf-input-output-unimplemented.md](../issues/M-io-vcf-input-output-unimplemented.md))
- [/] Auspice JSON (substitutions implemented; required `meta.updated` is absent and inference metadata remains incomplete: [kb/issues/H-io-auspice-v2-required-updated-missing.md](../issues/H-io-auspice-v2-required-updated-missing.md), [kb/issues/M-timetree-tree-output-inference-metadata-incomplete.md](../issues/M-timetree-tree-output-inference-metadata-incomplete.md))
- [ ] Skyline TSV/plot
- [ ] Substitution rates TSV
- [ ] Outliers TSV
- [ ] Branch mutations TXT
- [ ] Sequence evolution model TXT

## v1-Only Formats

- [x] PhyloXML
- [/] UShER MAT (partial - reference nucleotides can use the parent allele, and ancestral parsimony cannot supply TreeIR; [kb/issues/H-io-usher-ref-nuc-uses-parent-allele.md](../issues/H-io-usher-ref-nuc-uses-parent-allele.md), [kb/issues/N-ancestral-auspice-json-not-produced.md](../issues/N-ancestral-auspice-json-not-produced.md))
- [x] YAML serialization
- [x] Compressed FASTA output
- [x] Streaming readers/writers with automatic decompression
