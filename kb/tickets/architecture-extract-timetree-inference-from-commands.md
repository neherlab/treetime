# Extract timetree inference and optimization from commands/ to top-level module

## Description

Create top-level `src/timetree/` domain module. Move inference (479 lines), optimization (802 lines), and supporting files from `commands/timetree/`. After extraction, `commands/timetree/` drops from ~4900 to ~1900 lines.

## What to move

### inference/ (479 lines)

- `backward_pass.rs` (90 lines)
- `forward_pass.rs` (97 lines)
- `branch_length_likelihood.rs` (89 lines)
- `runner.rs` (195 lines) - `run_timetree()` core function
- `__tests__/` - inference tests

### optimization/ (802 lines)

- `polytomy.rs` (469 lines)
- `relaxed_clock.rs` (125 lines)
- `clock_filter.rs` (104 lines)
- `reroot.rs` (97 lines)
- `__tests__/` - optimization tests

### Supporting files

- `utils.rs` (118 lines) - `initialize_node_divergences`, consumed by inference
- `initialization.rs` (126 lines) - partition setup, no args coupling
- `refinement.rs` (105 lines) - iteration refinement helpers

## What stays in commands/timetree/

- `args.rs` (307 lines) - CLI argument structs
- `run.rs` (571 lines) - command orchestration
- `output/` (603 lines) - auspice JSON, confidence intervals, dates, plots
- `convergence/` (334 lines) - iteration metrics

## Related issues

- Source: [M-core-remaining-architectural-debt-after-extraction](../issues/M-core-remaining-architectural-debt-after-extraction.md)
