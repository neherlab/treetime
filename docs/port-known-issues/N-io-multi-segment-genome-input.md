# Multi-segment genome input not wired

The design document (`docs/algorithms/sequence_evolution.md:130-131`) asks how flu genome segments should be handled: "should these sequences be saved as list of multiple sequences, or concatenated?"

## Current state

The partition architecture supports multiple independent partitions on the same tree, each with its own alphabet and model. Discrete traits (mugration) demonstrate this capability. Multi-segment genome input (loading separate FASTA files per segment as separate partitions) has no CLI wiring.

## Workaround

Concatenate segment alignments into a single FASTA. This loses per-segment model assignment but preserves all sequence data.

## Design document context

`docs/algorithms/sequence_evolution.md:130-131`: "For flu, genomes come in segments -- should these sequences be saved as list of multiple sequences, or concatenated? Flu genome segments all use the same alphabet, so this could be implemented either way."
