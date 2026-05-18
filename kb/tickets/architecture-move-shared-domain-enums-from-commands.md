# Move shared domain enums from commands/shared/ to domain modules

## Description

`commands/shared/args.rs` contains three domain enums with `clap::ValueEnum` derives:

- `BranchLengthMode` (`Input`, `Marginal`): used by optimize and timetree commands
- `RerootMode` (`LeastSquares`, `MinDev`, `Oldest`, `ClockFilter`, `Mrca`): used by clock and timetree commands
- `MethodAncestral` (`Joint`, `Marginal`, `Parsimony`): used by ancestral and homoplasy commands

These are domain concepts (how to treat branch lengths, how to root the tree, which reconstruction method to use), not CLI concepts. They belong in domain modules.

### Target locations

- `BranchLengthMode` to `optimize/params.rs` (alongside `BranchOptMethod`, `InitialGuessMode`)
- `RerootMode` to `clock/find_best_root/params.rs` (alongside `BranchPointOptimizationParams`)
- `MethodAncestral` to `ancestral/` (new file or existing params module)

Strip `clap::ValueEnum` derives. Replace with `strum::EnumString` + `strum::Display` for string parsing. CLI args files use clap's `value_parser` or wrapper types.

Delete `commands/shared/args.rs` and `commands/shared/mod.rs` after all enums moved.

## Validation

- Domain modules compile without `commands/` imports
- All existing tests pass
- CLI enum parsing behavior unchanged

## Related issues

Source: [H-core-multi-client-architecture-library-purity](../issues/H-core-multi-client-architecture-library-purity.md)
