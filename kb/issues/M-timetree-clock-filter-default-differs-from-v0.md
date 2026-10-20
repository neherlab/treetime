# Timetree clock-filter default is 3.0; v0 uses 4.0

v0's `--clock-filter` defaults to 4.0 interquartile ranges ([packages/legacy/treetime/treetime/argument_parser.py#L154-L161](../../packages/legacy/treetime/treetime/argument_parser.py#L154-L161)). v1 sets the default to 3.0 in both commands that filter clock outliers:

- `timetree`: field `clock_filter` [packages/app-commands/src/commands/timetree/args.rs#L321-L328](../../packages/app-commands/src/commands/timetree/args.rs#L321-L328)
- `clock`: field `clock_filter` [packages/app-commands/src/commands/clock/args.rs#L125-L131](../../packages/app-commands/src/commands/clock/args.rs#L125-L131)

No decision records the change.

## Impact

With default flags, v1 marks more tips as clock outliers than v0, which changes the clock rate, the root and the node dates.

## Fix direction

Set the default to 4.0 in both argument structs, then regenerate the schemas and the CLI reference with `just gen`.

## Related issues

- [M-clock-filter-residual-parity.md](M-clock-filter-residual-parity.md): the residual computation itself differs. Both changes are necessary before the outlier sets can match v0
