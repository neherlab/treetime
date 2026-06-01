# Rename ClockParams\_ to avoid underscore suffix

## Description

`packages/treetime/src/clock/pipeline.rs` defines `pub struct ClockParams_` with an underscore suffix to avoid collision with `ClockParams` from `clock/clock_regression.rs`. The underscore naming is a hack. Rename to `ClockPipelineParams` or restructure imports.

## Related issues

- Source: [H-core-command-module-shared-ops-entanglement.md](../issues/H-core-command-module-shared-ops-entanglement.md)
