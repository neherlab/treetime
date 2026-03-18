# Optimize command accepts only a single alignment

The optimize command takes one `--aln` argument. The design document (`docs/algorithms/optimize.md:6,9`) specifies multiple alignments as partitions and a config file format for complex inputs with "multiple alignments and models along with discrete characters."

## Current state

`packages/treetime/src/commands/optimize/args.rs` accepts a single `--aln` path. The command builds one sparse and one dense partition from the same alignment. The partition architecture supports multiple independent partitions, but no CLI path loads multiple alignments with per-partition model assignment.

## Design document specification

`docs/algorithms/optimize.md:6`: "alignment(s) (corresponding to partitions)"
`docs/algorithms/optimize.md:9`: "For more complex inputs (multiple alignments and models along with discrete characters we might need a config file format)."

## Workaround

Concatenate alignments before running optimize (loses per-partition model control).
