# Extract coalescent/ from commands/timetree/ to top-level module

## Description

Move `commands/timetree/coalescent/` (1444 lines, 13 files) to top-level `src/coalescent/`. Self-contained domain logic with zero imports from `commands::timetree::args` or `commands::timetree::run`. Mechanical path rewriting, same pattern as ancestral/clock/optimize extraction.

## What to move

- `coalescent.rs` (79 lines) - coalescent model
- `skyline.rs` (315 lines) - skyline fitting
- `optimize_tc.rs` (161 lines) - Tc optimization
- `contributions.rs` (191 lines) - coalescent contributions
- `edge_data.rs` (146 lines) - edge data for coalescent
- `events.rs` (55 lines) - event collection
- `integration.rs` (94 lines) - numerical integration
- `lineage_dynamics.rs` (43 lines) - lineage tracking
- `piecewise_constant_fn.rs` (72 lines) - piecewise constant functions
- `piecewise_linear_fn.rs` (108 lines) - piecewise linear functions
- `time_coordinate.rs` (116 lines) - time coordinate utilities
- `total_lh.rs` (49 lines) - total likelihood
- `__tests__/` (1192 lines) - coalescent tests

## Related issues

- Source: [M-core-remaining-architectural-debt-after-extraction](../issues/M-core-remaining-architectural-debt-after-extraction.md)
