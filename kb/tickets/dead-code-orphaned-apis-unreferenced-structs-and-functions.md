# Remove orphaned APIs: unreferenced structs, functions, re-exports, and marker trait

## Description

Four instances of dead code that increase binary size and cognitive load.

## Instances

### MugrationParams and params() unreferenced

`packages/treetime/src/commands/mugration/input.rs:42-72:`

`struct MugrationParams` and its accessor `fn params()` are defined but never referenced anywhere in the workspace. Remove both.

### collect_edge_contributions -- pub but zero callers

`packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs:13:`

Public function with no callers in the workspace. Remove or mark as `#[allow(dead_code)]` with justification if needed for future use.

### lib.rs re-exports with no workspace consumer

`packages/treetime/src/lib.rs:15-17:`

Re-exports `treetime_graph`, `treetime_io`, `treetime_primitives` but no workspace crate imports them through this path. Remove the re-exports.

### PartitionMarginal -- empty marker trait

`packages/treetime/src/representation/partition/traits.rs:98:`

`trait PartitionMarginal` has no methods and no associated types. Used as a supertrait bound on `trait PartitionMarginalOps<N, E>`. The marker adds no contract beyond requiring `Send + Sync`, which could be expressed as direct bounds. Remove the marker trait and replace with direct `Send + Sync` bounds on `PartitionMarginalOps`.

## Related issues

- Source: [N-dead-code-orphaned-apis.md](../issues/N-dead-code-orphaned-apis.md) -- delete after full resolution
