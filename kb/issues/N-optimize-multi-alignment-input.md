# Optimize command loads all alignment files into one partition

> [!WARNING]
> **Needs review.** The command now accepts several alignment files: `pub alignment: Vec<PathBuf>` in `struct AlignmentArgs` [packages/app-commands/src/commands/shared/alignment.rs#L38](../../packages/app-commands/src/commands/shared/alignment.rs#L38). The records of all files form one alignment, so the files do not become separate partitions. The limitation below is restated for this behavior; confirm that the scope is still correct.

The optimize command reads every `--alignment` (alias `--aln`) file into one alignment and one partition. The intended input is one partition per alignment, each with its own substitution model, and a config file format for complex inputs that combine several alignments, models, and discrete characters ([kb/proposals/config-file-multi-partition.md](../proposals/config-file-multi-partition.md)).

## Current state

`struct AlignmentArgs` [packages/app-commands/src/commands/shared/alignment.rs#L21-L39](../../packages/app-commands/src/commands/shared/alignment.rs#L21-L39) accepts several FASTA files, and `fn read_alignment()` [packages/app-commands/src/commands/shared/alignment.rs#L52-L58](../../packages/app-commands/src/commands/shared/alignment.rs#L52-L58) appends their records into one list. The optimize command [packages/app-commands/src/commands/optimize/run.rs#L41](../../packages/app-commands/src/commands/optimize/run.rs#L41) builds one partition from it, dense or sparse as selected by `fn Representation::resolve()` [packages/treetime/src/partition/create.rs#L25-L31](../../packages/treetime/src/partition/create.rs#L25-L31). The partition architecture supports multiple independent partitions, but no CLI path loads multiple alignments with per-partition model assignment.

## Workaround

Concatenate alignments column-wise before running optimize (loses per-partition model control).

## Related

### Known issues

- [M-io-sequence-name-matching-unreliable.md](M-io-sequence-name-matching-unreliable.md) -- name matching affects multi-alignment attachment
- [N-io-multi-segment-genome-input.md](N-io-multi-segment-genome-input.md) -- related multi-input limitation

### Proposals

- [kb/proposals/config-file-multi-partition.md](../proposals/config-file-multi-partition.md) -- configuration file format for multiple alignments
- [kb/proposals/unified-input-format-support.md](../proposals/unified-input-format-support.md) -- alternative input via formats with embedded sequences
