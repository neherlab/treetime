# Timetree clock-filter default is 3.0; v0 uses 4.0

v0's `--clock-filter` defaults to 4.0 interquartile ranges ([packages/legacy/treetime/treetime/argument_parser.py#L155-L161](../../packages/legacy/treetime/treetime/argument_parser.py#L155-L161)). v1's timetree default is 3.0 ([packages/app-cli/src/commands/timetree/args.rs](../../packages/app-cli/src/commands/timetree/args.rs), field `clock_filter`), and the `clock` command also uses 3.0. No decision records the change.

## Impact

With default flags, v1 marks more tips as clock outliers than v0, which changes the clock rate, the root and the node dates.

## Related issues

- [H-cli-timetree-config-disables-clock-filter.md](H-cli-timetree-config-disables-clock-filter.md): config files set the value to 0.0
- [M-clock-filter-residual-parity.md](M-clock-filter-residual-parity.md): the residual computation itself differs
