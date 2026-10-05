# Multi-segment genome input not wired

Segmented genomes such as influenza need a representation: one partition per segment, or one concatenated sequence. All flu segments use the same alphabet, so either works for sequences. Discrete traits such as location or host have their own alphabets, which favors the partition structure.

## Current state

The partition architecture supports multiple independent partitions on the same tree, each with its own alphabet and model. Discrete traits (mugration) demonstrate this capability. Multi-segment genome input (loading separate FASTA files per segment as separate partitions) has no CLI wiring.

The public alignment-input design discussion proposes separate nucleotide and amino-acid alignment partitions [[issue](https://github.com/neherlab/treetime/issues/306)]. A maintainer response identifies sequence alphabets, discrete traits, sparse/dense storage, and configuration-based partition specifications as parts of the same partition abstraction [[comment](https://github.com/neherlab/treetime/issues/306#issuecomment-2565383822)].

The branch `worktree/feat/multi-segment-genome-input` has 9 commits implementing a `--segment` flag, but this work has not been merged.

## v0 comparison

v0 accepts only a single alignment file. Multi-segment genomes are handled by concatenation before running TreeTime. The `--aln` flag accepts one file.

## Workaround

Concatenate segment alignments into a single FASTA. This loses per-segment model assignment but preserves all sequence data.

## Related

### Known issues

- [M-io-sequence-name-matching-unreliable](M-io-sequence-name-matching-unreliable.md) -- name matching affects multi-segment attachment
- [N-io-large-dataset-memory-constraint](N-io-large-dataset-memory-constraint.md) -- memory constraints compound with multiple segments

### Proposals

- [config-file-multi-partition](../proposals/config-file-multi-partition.md) -- configuration file format that would subsume per-segment CLI flags
- [unified-input-format-support](../proposals/unified-input-format-support.md) -- alternative input via formats with embedded sequences
