# Dead code and orphaned APIs

## Summary

Five instances of dead code: unreferenced structs, pub functions with zero callers, re-exports with no consumers, and a marker trait with no contract.

## Instances

### MugrationParams and MugrationInput::params() unreferenced

`packages/treetime/src/commands/mugration/input.rs:42-72:`

`struct MugrationParams` and its accessor `fn params()` are defined but never referenced anywhere in the workspace.

### collect_edge_contributions pub but zero callers

`packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs:13:`

Public function with no callers in the workspace.

### lib.rs pub use re-exports with no workspace consumer

`packages/treetime/src/lib.rs:15-17:`

Re-exports `treetime_graph`, `treetime_io`, `treetime_primitives` but no workspace crate imports them through this path.

### DenseSeqDis::new() hard-codes log_lh: 0.0

`packages/treetime/src/representation/payload/dense.rs:53-55:`

`struct DenseSeqDis` holds a 2D probability array and a `log_lh: f64` field tracking the log-likelihood normalization constant. The constructor sets `log_lh: 0.0` with no parameter to override. Callers that need a non-zero initial `log_lh` must mutate the field after construction, violating the "initialize fully upfront" principle.

### PartitionMarginal empty marker trait

`packages/treetime/src/representation/partition/traits.rs:98:`

`trait PartitionMarginal` has no methods and no associated types. Used as a supertrait bound on `trait PartitionMarginalOps<N, E>`. The marker adds no contract beyond requiring `Send + Sync`, which could be expressed as direct bounds.

## Impact

- Dead code increases binary size and cognitive load
- Marker trait without contract creates false abstraction boundary
- Constructor requiring post-construction mutation is error-prone
