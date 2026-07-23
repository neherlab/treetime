# treetime crate README describes a stale module layout

[`packages/treetime/README.md`](../../packages/treetime/README.md) still lists a nonexistent `cli` module, describes `commands` as analysis pipelines, and names the partition and payload layer `representation` [`packages/treetime/README.md#L17-L31`](../../packages/treetime/README.md#L17-L31).

## Source mismatch

- CLI entry points live in `app-cli`, not a `treetime::cli` module.
- Scientific pipelines live in domain modules such as ancestral, optimize, and timetree; `commands` owns argument and I/O orchestration around them.
- Partition storage/capabilities and graph payloads are separate `partition` and `payload` modules rather than one `representation` layer.
- CLI, HTTP, N-API, desktop, and web consumers make command ownership a multi-client question.

The README is an entry point for architectural navigation. Incorrect ownership descriptions send contributors toward nonexistent paths and reinforce an unapproved CLI-only boundary.

## Required update

Correct the module inventory for paths and responsibilities that are already factual. Describe the application boundary as unresolved and link the governing issue rather than selecting one ownership option in documentation.

## Validation

- Every named module resolves in `packages/treetime/src/lib.rs` or is clearly identified as another workspace crate.
- Descriptions match imports and public exports in source.
- The README contains no proposed ownership presented as current architecture.

## Related issues

- [H-core-command-module-shared-ops-entanglement.md](H-core-command-module-shared-ops-entanglement.md)
- [H-core-multi-client-architecture-library-purity.md](H-core-multi-client-architecture-library-purity.md)
