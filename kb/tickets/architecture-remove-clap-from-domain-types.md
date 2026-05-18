# Remove clap derives from domain types

## Description

6 domain files derive `clap::ValueEnum` or `clap::Args`. This prevents the `treetime` library crate from compiling without the `clap` feature. Currently the feature gate is cosmetic: `treetime-cli` always enables `features = ["clap"]`, so the separation is never tested.

### Domain files with clap derives

- `alphabet/alphabet.rs`: `AlphabetName` derives `ValueEnum`
- `clock/clock_regression.rs`: type with clap derive (verify exact type)
- `clock/find_best_root/params.rs`: `BranchPointOptimizationParams` or related type
- `gtr/get_gtr.rs`: `GtrModelName` derives `ValueEnum`
- `optimize/params.rs`: `BranchOptMethod`, `InitialGuessMode` derive `ValueEnum`
- `seq/gap_fill.rs`: `GapFillMode` derives `ValueEnum`

### commands/shared/ enums

`commands/shared/args.rs` contains `BranchLengthMode`, `RerootMode`, `MethodAncestral` with `ValueEnum` derives. These are domain enums that should move to their respective domain modules before this ticket.

### Approach

Replace `clap::ValueEnum` with `strum::EnumString` + `strum::Display` (or `strum::AsRefStr`) for string parsing. In CLI args files, use clap's `value_parser` with the strum-derived `FromStr` impl, or create thin wrapper enums with `ValueEnum`.

### Prerequisite

`commands/shared/args.rs` domain enums moved to domain modules (separate ticket).

## Validation

- `treetime` crate compiles without `clap` feature: `cargo check -p treetime --no-default-features`
- `treetime-cli` builds and all tests pass
- CLI enum parsing behavior unchanged (same accepted values, same error messages)

## Related issues

Source: [H-core-multi-client-architecture-library-purity](../issues/H-core-multi-client-architecture-library-purity.md)
