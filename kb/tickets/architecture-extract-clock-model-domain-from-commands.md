# Extract clock model domain from commands/clock/ to domain module

## Description

`commands/clock/` contains 11 non-command files (1443 lines total) defining the clock model domain. These are consumed by `timetree` (20 import sites), `representation/payload/`, `representation/partition/`, and `cli/`. The clock model is domain logic, not command logic.

## What to move

11 files in `commands/clock/`:

- `clock_model.rs` -- `ClockModel`
- `clock_set.rs` -- `ClockSet`
- `clock_traits.rs` -- `ClockNode`, `ClockEdge`
- `clock_regression.rs` -- `ClockParams`, `estimate_clock_model_with_reroot()`, `estimate_clock_model_with_reroot_policy()`
- `clock_filter.rs` -- `clock_filter_inplace()`
- `clock_output.rs` -- `write_clock_model()`
- `reroot.rs` -- `RerootChanges`, `RerootParams`
- `rtt.rs` -- root-to-tip distance
- `date_constraints.rs` -- `load_date_constraints()`
- `clock_graph.rs` -- clock graph types
- `assign_dates.rs` -- date assignment

Total: 1443 lines.

## Consumers

### timetree (20 import sites)

`ClockModel`, `ClockSet`, `ClockNode`, `ClockEdge`, `ClockParams`, `RerootChanges`, `RerootParams`, `BranchPointOptimizationParams`, `clock_filter_inplace()`, `estimate_clock_model_with_reroot()`, `estimate_clock_model_with_reroot_policy()`, `load_date_constraints()`, `write_clock_model()`

### representation/ (reverse dependency -- core importing from commands)

- `representation/payload/timetree.rs#L1-L4` -> `ClockSet`, `ClockEdge`, `ClockNode`, `DateConstraintNode`, `TimetreeEdge`, `TimetreeNode` from `commands/clock/` and `commands/timetree/`

### representation/partition/ (reverse dependency)

- `representation/partition/timetree.rs#L1` -> `PartitionRerootOps`, `PartitionTimetreeAll` from `commands/timetree/partition_ops`

### cli/ (reverse dependency)

- `cli/rtt_chart.rs#L2-L3`, `cli/rtt_chart_render.rs#L1-L2` -> `ClockModel`, `ClockRegressionResult` from `commands/clock/`

## Target location

Extract to a new domain module outside `commands/` (e.g. `src/clock/` or a new crate). The `commands/clock/` command module should become a thin wrapper.

## Related issues

- Source: [H-core-command-module-shared-ops-entanglement.md](../issues/H-core-command-module-shared-ops-entanglement.md) -- delete after full resolution
