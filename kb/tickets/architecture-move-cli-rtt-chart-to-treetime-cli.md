# Move cli/rtt_chart modules to treetime-cli

## Description

`cli/rtt_chart.rs` (177 lines) and `cli/rtt_chart_render.rs` (174 lines) render root-to-tip regression charts. Called exclusively by `commands/clock/run.rs`. They import from `clock/` domain module, not from commands. Move them to `treetime-cli` alongside the commands that use them.

## What to change

- Move `packages/treetime/src/cli/rtt_chart.rs` to `packages/treetime-cli/src/cli/rtt_chart.rs`
- Move `packages/treetime/src/cli/rtt_chart_render.rs` to `packages/treetime-cli/src/cli/rtt_chart_render.rs`
- Update `commands/clock/run.rs` import from `crate::cli::` to the new location
- Remove `cli/` module from `treetime` crate if empty after move (only `mod.rs` remains)

## Related issues

- Source: [M-core-cli-library-separation-blockers](../issues/M-core-cli-library-separation-blockers.md) (B2)
