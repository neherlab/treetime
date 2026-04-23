# RTT chart aborts clock command on non-TTY output

`print_clock_regression_chart` at [packages/treetime/src/cli/rtt_chart.rs#L91](../../packages/treetime/src/cli/rtt_chart.rs#L91) calls `terminal::size()?` to determine the terminal width for chart rendering. When stdout is not connected to a TTY (piped output, CI environments, redirected to file), `terminal::size()` returns an IO error. The `?` operator propagates this error, aborting the entire `clock` command.

## Current mitigation

The call site at [packages/treetime/src/commands/clock/run.rs#L106-L108](../../packages/treetime/src/commands/clock/run.rs#L106-L108) is now gated on `is_tty()`, so the function is never called in non-TTY contexts. The abort described above does not occur in practice. The function itself still uses `terminal::size()?` without a fallback.

## Remaining fix

Replace `terminal::size()?` with `terminal::size().unwrap_or((80, 24))` inside `print_clock_regression_chart` to make the function itself resilient to non-TTY contexts, independent of call-site guards.
