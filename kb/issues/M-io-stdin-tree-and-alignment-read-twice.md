# Tree and alignment can both read standard input

A command accepts `-` as `--tree` and as an `--alignment` file in the same run. Both inputs then read standard input, and the run fails with an error that names the wrong cause.

## Current state

- `fn AlignmentArgs::paths()` rejects `-` when the `--alignment` list holds it more than once [packages/app-commands/src/commands/shared/alignment.rs#L53-L59](../../packages/app-commands/src/commands/shared/alignment.rs#L53-L59), but no check compares the alignment paths with `--tree`
- `fn read_input_tree()` reads the tree through `fn tree_read_file()`, which reads `-` from standard input to its end [packages/app-commands/src/commands/shared/tree_input.rs#L25-L26](../../packages/app-commands/src/commands/shared/tree_input.rs#L25-L26) [packages/treetime-io/src/tree.rs#L11-L17](../../packages/treetime-io/src/tree.rs#L11-L17)
- Every command that reads a tree and an alignment reads the tree first, so the alignment reads an empty standard input

## Reproduction

With the zika/20 dataset (`data/zika/20/`, alignment decompressed to `aln.fasta`):

- `treetime optimize --tree - --alignment - --output-all out < data/zika/20/tree.nwk` stops with "More than one third of terminal nodes (20 of 20) cannot be assigned a sequence. Are you sure the alignment belongs to the tree?"
- `treetime optimize --tree - --alignment - --output-all out < aln.fasta` stops with a Newick parse error at line 2, column 1
- Concatenating the tree and the FASTA into standard input stops with "The input contains more than one tree"

## Fix direction

Reject `-` in more than one input of a run before any input is read, with one message that names the inputs, in the place that owns the input paths of a command. The same rule then covers other inputs that accept `-`.

> [!IMPORTANT]
> **Investigation required.** The reproduction covers `--tree` and `--alignment` of `optimize`. Check which other inputs accept `-` (for example the metadata and dates tables of `timetree`, `clock` and `mugration`) and whether they can collide in the same way.

## Related

- [kb/features/io.md](../features/io.md) -- standard input entry
