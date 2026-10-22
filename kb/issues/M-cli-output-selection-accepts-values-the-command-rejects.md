# `--output-selection` accepts values that the command rejects at runtime

The `command_outputs!` macro adds the full tree-format set (`nwk`, `nexus`, `auspice`, `mat-pb`, `mat-json`, `graph-json`, `dot`) to the selection enum of every command ([packages/app-commands/src/commands/shared/output_args.rs#L27-L82](../../packages/app-commands/src/commands/shared/output_args.rs#L27-L82)). Which outputs a command can write is decided separately, by the output table of the command (`CommandKind::outputs()` in `packages/app-output/src/command_outputs.rs`) and by the output conditions of a run (`OutputCondition`, `UnavailableOutput` in `packages/app-output/src/output_conditions.rs`). A value that the parser accepts can therefore fail only after the inputs are read. The per-file output flags have the same drift when their writer is absent from the output table of the command.

> [!IMPORTANT]
> **Investigation required.** Every command lists the full tree-format set in its output table today, so the values that fail at runtime come from the output conditions. List each command and selection value that the parser accepts and the run rejects, with the error it prints.

## Fix direction

Generate the selection enum and the visible per-file flags of each command from the outputs that the command can write, and add tests that the parser rejects the other values.

## Affected code

- `packages/app-commands/src/commands/shared/output_args.rs`: `command_outputs!` and the selection enums
- `packages/app-output/src/command_outputs.rs`: the output table of each command
- `packages/app-output/src/output_conditions.rs`: the conditions that reject an output at runtime
