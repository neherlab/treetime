# RTT chart aborts clock command on non-TTY output

`print_clock_regression_chart` at [packages/treetime/src/cli/rtt_chart.rs#L89-L91](../../packages/treetime/src/cli/rtt_chart.rs#L89-L91) calls `terminal::size()?` to determine the terminal width for chart rendering. When stdout is not connected to a TTY (piped output, CI environments, redirected to file), `terminal::size()` returns an IO error. The `?` operator propagates this error, aborting the entire `clock` command.

## Trigger conditions

The function is called from [packages/treetime/src/commands/clock/run.rs#L103-L107](../../packages/treetime/src/commands/clock/run.rs#L103-L107), gated on `args.plot`. Running the clock command with `--plot-rtt` in a non-TTY context (e.g., `treetime clock ... | tee output.log`) causes the command to fail with an IO error.

## Fix

Replace `terminal::size()?` with `terminal::size().unwrap_or((80, 24))` to fall back to a default terminal width, or check `crossterm::tty::IsTty` before attempting chart rendering and skip the chart with a warning when output is not a TTY.
