# Move commands/ orchestration shell from library to CLI crate

## Description

After domain logic extraction (ancestral, clock, optimize, coalescent, timetree inference, mugration), the `commands/` module in the `treetime` library contains only thin CLI orchestration: argument structs (`args.rs`), command runners (`run.rs`), and output formatting (`output/`). Move this entire module to `treetime-cli`.

## What to move

All of `packages/treetime/src/commands/` to `packages/treetime-cli/src/commands/`:

- `ancestral/` - args.rs, run.rs
- `clock/` - args.rs, run.rs
- `homoplasy/` - args.rs, run.rs
- `mugration/` - args.rs, run.rs
- `optimize/` - args.rs, run.rs
- `prune/` - args.rs, run.rs
- `timetree/` - args.rs, run.rs, output/, initialization.rs, refinement.rs
- `mod.rs`

## What stays in treetime library

Nothing from `commands/`. The library crate exposes domain modules only:

- `alphabet/`, `ancestral/`, `clock/`, `coalescent/`, `constants/`, `gtr/`, `mugration/`, `optimize/`, `partition/`, `payload/`, `prune/`, `seq/`, `timetree/`

## Prerequisites

All of these must be completed first:

- `architecture-extract-timetree-inference-from-commands.md`
- ~~`architecture-extract-mugration-domain-logic.md`~~ (resolved)
- `architecture-break-representation-gtr-cycle.md`
- `architecture-break-representation-ancestral-cycle.md`

## Validation

- `treetime` library crate has no `commands/` module
- `treetime-cli` builds and passes all existing tests
- No circular dependencies between `treetime` and `treetime-cli`

## Related issues

Source: [kb/issues/H-core-multi-client-architecture-library-purity.md](../issues/H-core-multi-client-architecture-library-purity.md)
