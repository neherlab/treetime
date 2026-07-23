# Command output ownership is scattered

The shared output module owns selection, default sets, path resolution, and tree-format dispatch [`packages/treetime/src/commands/shared/output.rs#L28`](../../packages/treetime/src/commands/shared/output.rs#L28) [`packages/treetime/src/commands/shared/tree_output.rs#L72`](../../packages/treetime/src/commands/shared/tree_output.rs#L72). Each command runner still registers non-tree fields, checks data prerequisites, chooses projection inputs and comment providers, and writes non-tree outputs independently.

`OutputCoreArgs::resolve()` is pure, and tree projection and writing are centralized. Ownership of a complete output plan and materialization contract remains distributed across command modules.

## Impact

Adding an output kind or changing its data prerequisites requires inspecting the selection registry, command capability/default sets, one or more command runners, and format-specific writers.

Each runner supplies another keyed slice of non-tree fields and chooses when required domain facts are available. Tree output is centralized, but non-tree output still requires command-specific checks and dispatch. The effective output contract is therefore the intersection of registry declarations and runner control flow.

## Design question

Define whether output descriptors own command capability and data prerequisites, and whether analysis results expose renderable in-memory views that a shared materializer consumes. Filesystem adapters must remain separable from pure planning and projection tests.

## Required properties

- One descriptor states selection name, extension, command availability, default status, and required domain facts.
- Resolution produces a pure complete plan without creating directories or files.
- Materialization receives analysis facts and chooses projection and writer behavior in one boundary.
- In-memory projections remain testable without filesystem I/O.
- Missing prerequisites fail before partial output is written.

## Validation

- Registry tests prove every selection has one complete descriptor and every command default is supported.
- Planning tests cover explicit paths, output directories, compression, conflicts, and unavailable outputs without touching the filesystem.
- Materialization tests cover each output kind from typed analysis results.
- Integration tests verify transactional behavior when one selected writer fails.

## Related issues

- [H-core-command-module-shared-ops-entanglement.md](H-core-command-module-shared-ops-entanglement.md)
