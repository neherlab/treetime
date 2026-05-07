# Move partition and payload traits from commands/timetree/ to representation/

## Description

Trait definitions in `commands/timetree/` are consumed by `representation/partition/` and `representation/payload/`, creating a reverse dependency where the representation layer imports from a command module. These traits belong in `representation/` alongside the types that implement them.

## What to move

- `packages/treetime/src/commands/timetree/partition_ops.rs` (60 lines) -- defines `PartitionRerootOps`, `PartitionTimetreeAll`
- `packages/treetime/src/commands/timetree/timetree_traits.rs` (18 lines) -- defines timetree traits

Total: 78 lines.

## Consumers (reverse dependencies)

- `representation/partition/marginal_dense.rs#L3` -> `PartitionRerootOps`, `PartitionTimetreeAll` from `commands/timetree/partition_ops`
- `representation/partition/marginal_sparse.rs#L4` -> `PartitionRerootOps`, `PartitionTimetreeAll` from `commands/timetree/partition_ops`
- `representation/partition/timetree.rs#L1` -> `PartitionRerootOps`, `PartitionTimetreeAll` from `commands/timetree/partition_ops`
- `representation/payload/timetree.rs#L1-L4` -> `TimetreeEdge`, `TimetreeNode` from `commands/timetree/`

## Target location

Move to `representation/` where the implementing types live (e.g. `representation/partition/` or `representation/traits/`).

## Related issues

- Source: [H-core-command-module-shared-ops-entanglement.md](../issues/H-core-command-module-shared-ops-entanglement.md) -- delete after full resolution
