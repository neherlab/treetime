# Tree and alignment can both read standard input

A command accepts `-` as `--tree` and as an `--alignment` file in the same run. Both inputs then read standard input, and the run fails with an error that names the wrong cause.

## Current state

- `fn AlignmentArgs::paths()` rejects `-` when the `--alignment` list holds it more than once [packages/app-commands/src/commands/shared/alignment.rs#L55-L61](../../packages/app-commands/src/commands/shared/alignment.rs#L55-L61), but no check compares the alignment paths with `--tree`
- `fn read_input_tree()` reads the tree through `fn tree_read_file()`, which reads `-` from standard input to its end [packages/app-commands/src/commands/shared/tree_input.rs#L25-L26](../../packages/app-commands/src/commands/shared/tree_input.rs#L25-L26) [packages/treetime-io/src/tree.rs#L11-L17](../../packages/treetime-io/src/tree.rs#L11-L17)
- Every command that reads a tree and an alignment reads the tree first, so the alignment reads an empty standard input

## Reproduction

With the zika/20 dataset (`data/zika/20/`, alignment decompressed to `aln.fasta`):

- `treetime optimize --tree - --alignment - --output-all out < data/zika/20/tree.nwk` stops with "More than one third of terminal nodes (20 of 20) cannot be assigned a sequence. Are you sure the alignment belongs to the tree?"
- `treetime optimize --tree - --alignment - --output-all out < aln.fasta` stops with a Newick parse error at line 2, column 1
- Concatenating the tree and the FASTA into standard input stops with "The input contains more than one tree"

## Fix direction

Reject `-` in more than one input of a run before any input is read, with one message that names the inputs, in the place that owns the input paths of a command. The same rule then covers other inputs that accept `-`.

## Inputs that read `-` from standard input

Every input read through `fn read_file_with()` of `treetime-utils` reads `-` from standard input, and only the `--alignment` list checks for a second `-`:

- `--tree` of every command (`fn tree_read_file()`)
- `--alignment` of `ancestral`, `optimize`, `timetree`, `prune`, `homoplasy`
- `--metadata` of `clock`, `timetree`, `mugration`, and `--weights` of `mugration`
- `--drms` of `homoplasy`
- `--prune-nodes-list-file` of `prune`
- `--topology-order-target-file` of `optimize`, `timetree`, `prune`, `mugration`, `homoplasy`
- `--config` of every command

Each pair of these inputs in one run can collide in the same way. `--annotation` and `--aa-root-sequence` of `ancestral` do not take part, because their existence check rejects `-` before the read [packages/app-commands/src/commands/ancestral/aa_node_data.rs#L47-L48](../../packages/app-commands/src/commands/ancestral/aa_node_data.rs#L47-L48).

## Related

- [kb/features/io.md](../features/io.md) -- standard input entry
