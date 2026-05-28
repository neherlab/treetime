# Extract shared gap-fill CLI argument struct

Three command argument structs implement identical `effective_gap_fill()` methods.

## Current state

`commands/ancestral/args.rs:133-139`, `commands/optimize/args.rs:127-133`, and `commands/timetree/args.rs:310-316` each define `effective_gap_fill()` with the same `if keep_overhangs { None } else { gap_fill }` logic.

## Target state

A shared `GapFillArgs { gap_fill: GapFill, keep_overhangs: bool }` struct in `commands/shared/` with `effective_gap_fill()` implemented once. Each command embeds via `#[command(flatten)]`.

## Implementation

1. Define `GapFillArgs` in `commands/shared/args.rs` (or a new file) with clap derive attributes matching the current fields
2. Implement `effective_gap_fill()` on `GapFillArgs`
3. Replace the `gap_fill` and `keep_overhangs` fields in all three command arg structs with `#[command(flatten)] gap_fill_args: GapFillArgs`
4. Update all call sites from `args.effective_gap_fill()` to `args.gap_fill_args.effective_gap_fill()`
5. Verify CLI help text unchanged

## Related issues

Source: `kb/issues/N-core-gap-fill-args-duplication.md`
