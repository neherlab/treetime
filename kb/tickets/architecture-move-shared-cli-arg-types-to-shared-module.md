# Move shared CLI arg types to shared module

## Description

Three CLI arg types defined in individual command modules are consumed by other commands, creating cross-command dependencies. They should live in a shared location.

## What to move

- `BranchLengthMode` -- defined in `packages/treetime/src/commands/timetree/args.rs#L21`
- `RerootMode` -- defined in `packages/treetime/src/commands/timetree/args.rs#L40`
- `MethodAncestral` -- defined in `packages/treetime/src/commands/ancestral/args.rs#L11`

## Cross-command imports

- `commands/clock/args.rs` imports `BranchLengthMode` and `RerootMode` from `commands/timetree/args.rs`
- `commands/clock/args.rs` imports `MethodAncestral` from `commands/ancestral/args.rs`
- `commands/timetree/args.rs` imports `MethodAncestral` from `commands/ancestral/args.rs`

## Target location

Move to a shared args module (e.g. `commands/shared/args.rs` or `treetime-cli` crate) so that no command depends on another command's arg definitions.

## Related issues

- Source: [H-core-command-module-shared-ops-entanglement.md](../issues/H-core-command-module-shared-ops-entanglement.md) -- delete after full resolution
