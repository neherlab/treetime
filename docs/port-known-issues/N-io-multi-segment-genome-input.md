# Multi-segment genome input not wired

The design document ([docs/algorithms/sequence_evolution.md](../algorithms/sequence_evolution.md)) asks how flu genome segments should be handled: "should these sequences be saved as list of multiple sequences, or concatenated?"

## Current state

The partition architecture supports multiple independent partitions on the same tree, each with its own alphabet and model. Discrete traits (mugration) demonstrate this capability. Multi-segment genome input (loading separate FASTA files per segment as separate partitions) has no CLI wiring.

The branch `worktree/feat/multi-segment-genome-input` has 9 commits implementing a `--segment` flag, but this work has not been merged.

## v0 comparison

v0 accepts only a single alignment file. Multi-segment genomes are handled by concatenation before running TreeTime. The `--aln` flag accepts one file.

## Workaround

Concatenate segment alignments into a single FASTA. This loses per-segment model assignment but preserves all sequence data.

## Design document context

[docs/algorithms/sequence_evolution.md](../algorithms/sequence_evolution.md): "For flu, genomes come in segments -- should these sequences be saved as list of multiple sequences, or concatenated? Flu genome segments all use the same alphabet, so this could be implemented either way."

## Related

- [config-file-multi-partition](../port-proposals/config-file-multi-partition.md) -- proposal for configuration file format that would subsume per-segment CLI flags
- [docs/algorithms/sequence_evolution.md](../algorithms/sequence_evolution.md) -- source of the multi-segment requirement
