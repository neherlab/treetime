# Gap-fill argument logic duplicated across three commands

Three command argument structs implement identical `effective_gap_fill()` methods that return `GapFill::None` when `keep_overhangs` is set, otherwise return the configured `gap_fill` value.

## Locations

- `packages/treetime/src/commands/ancestral/args.rs:133-139:` `TreetimeAncestralArgs::effective_gap_fill()`
- `packages/treetime/src/commands/optimize/args.rs:127-133:` `TreetimeOptimizeArgs::effective_gap_fill()`
- `packages/treetime/src/commands/timetree/args.rs:310-316:` `TreetimeTimetreeArgs::effective_gap_fill()`

All identical:

```rust
fn effective_gap_fill(&self) -> GapFill {
  if self.keep_overhangs { GapFill::None } else { self.gap_fill }
}
```

## Impact

Trivial duplication. User-facing CLI policy logic that must stay consistent across commands.

## Action

Extract a shared clap-flattenable struct `GapFillArgs { gap_fill: GapFill, keep_overhangs: bool }` with `effective_gap_fill()` implemented once. Each command embeds the shared struct via `#[command(flatten)]`.
