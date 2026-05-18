# Remove clap derives from domain types

## Description

6 domain-layer files depend on `clap` for `ValueEnum` or `Args` derives. Domain types should not know about CLI parsing. Library consumers pull in `clap` as transitive dependency.

## Affected files

- `clock/clock_regression.rs` - `ClockParams` has `#[derive(Args)]`
- `clock/find_best_root/params.rs` - 5 types with `ValueEnum`/`Args`
- `optimize/args.rs` - `BranchOptMethod`, `InitialGuessMode` with `ValueEnum`
- `alphabet/alphabet.rs` - `AlphabetName` with `ValueEnum`
- `gtr/get_gtr.rs` - `GtrModelName` with `ValueEnum`
- `seq/gap_fill.rs` - `GapFillMode` with `ValueEnum`

## Fix

Create CLI wrapper types in `commands/` args files with `ValueEnum`/`Args` derives. Use `From` impls to convert to domain types. Domain types keep `Serialize`/`Deserialize` but drop `clap` dependency.

## Validation

- `clap` removed from `treetime` crate `[dependencies]` in `Cargo.toml`
- `cargo build` succeeds
- All existing tests pass
- CLI behavior unchanged

## Related issues

- Source: [M-core-remaining-architectural-debt-after-extraction](../issues/M-core-remaining-architectural-debt-after-extraction.md)
- Source: [H-core-multi-client-architecture-library-purity](../issues/H-core-multi-client-architecture-library-purity.md)
