# Optimize command accepts only a single alignment

The optimize command takes one `--aln` argument. The design document (`do../_raw/optimize.md:6,9`) specifies multiple alignments as partitions and a config file format for complex inputs with "multiple alignments and models along with discrete characters."

## Current state

`packages/treetime/src/commands/optimize/args.rs` accepts a single `--aln` path. The command builds one sparse and one dense partition from the same alignment. The partition architecture supports multiple independent partitions, but no CLI path loads multiple alignments with per-partition model assignment.

## Design document specification

`do../_raw/optimize.md:6`: "alignment(s) (corresponding to partitions)"
`do../_raw/optimize.md:9`: "For more complex inputs (multiple alignments and models along with discrete characters we might need a config file format)."

## Workaround

Concatenate alignments before running optimize (loses per-partition model control).

## Related

### Known issues

- [M-io-sequence-name-matching-unreliable](M-io-sequence-name-matching-unreliable.md) -- name matching affects multi-alignment attachment
- [M-io-sequence-attachment-quadratic](M-io-sequence-attachment-quadratic.md) -- attachment performance with multiple alignments
- [N-io-multi-segment-genome-input](N-io-multi-segment-genome-input.md) -- related multi-input limitation

### Proposals

- [config-file-multi-partition](../proposals/config-file-multi-partition.md) -- configuration file format for multiple alignments
- [unified-input-format-support](../proposals/unified-input-format-support.md) -- alternative input via formats with embedded sequences
