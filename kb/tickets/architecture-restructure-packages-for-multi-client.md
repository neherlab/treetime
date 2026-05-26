# Restructure existing packages for multi-client architecture

## Description

Reorganize the existing workspace so that the `treetime` library crate is a pure domain library consumable by any client (CLI, Python, desktop, WASM). The `treetime-cli` crate becomes one of several thin client crates.

This is an umbrella ticket coordinating prerequisite work items. Each item has its own ticket.

## Execution order

The work items have a dependency chain. Execute in this order:

### ~~1. Remove clap from library domain types~~ Cancelled

Decision: keep `clap` as a direct dependency of the `treetime` library crate.

### 2. Extract remaining domain logic from commands/

Three extraction tickets, independent of each other:

- ~~`architecture-extract-timetree-convergence-confidence.md`~~ (resolved: convergence and confidence extracted to `src/timetree/`)
- ~~`architecture-extract-mugration-domain-logic.md`~~ (resolved: mugration domain logic extracted to `src/mugration/`)
- ~~`architecture-extract-prune-domain-logic.md`~~ (resolved: pruning algorithms extracted to `src/prune/`)

### ~~3. Move shared domain enums~~ Done

- ~~`architecture-move-shared-domain-enums-from-commands.md`~~ (resolved: `BranchLengthMode`, `RerootMode`, `MethodAncestral` moved to domain modules)

### 4. Move commands/ shell to CLI crate

After steps 1-3, `commands/` contains only thin orchestration (`args.rs`, `run.rs`, `output/`). Move the entire `commands/` module to `treetime-cli`. The library crate no longer has a `commands/` module.

Ticket: `architecture-move-commands-to-cli-crate.md`

### 5. Add client crates

Two tickets, depend on steps 1-4:

- `architecture-add-treetime-python-crate.md` - PyO3 bindings + Python CLI
- `architecture-add-treetime-desktop-crate.md` - thin desktop wrapper

## Target architecture

```
treetime (library)
  +-- alphabet/
  +-- ancestral/
  +-- clock/
  +-- coalescent/
  +-- constants/
  +-- gtr/
  +-- mugration/
  +-- optimize/
  +-- partition/
  +-- seq/
  +-- timetree/

treetime-cli (client)
  +-- commands/        (moved from treetime, thin orchestration only)
  |   +-- ancestral/
  |   +-- clock/
  |   +-- homoplasy/
  |   +-- mugration/
  |   +-- optimize/
  |   +-- prune/
  |   +-- timetree/
  +-- bin/

treetime-python (client)
  +-- src/             (PyO3 bindings)
  +-- python/          (Python package + CLI)

treetime-desktop (client)
  +-- src/             (IPC handlers)
  +-- frontend/        (minimal UI, deferred)
```

## Validation

- `treetime` crate has no `commands/` module
- `treetime-cli` builds and passes all existing tests
- `treetime-python` builds a wheel and runs a smoke test
- `treetime-desktop` compiles (frontend deferred)

## Related issues

Source: [H-core-multi-client-architecture-library-purity](../issues/H-core-multi-client-architecture-library-purity.md)

Coordinates:

- ~~architecture-extract-timetree-convergence-confidence.md~~ (resolved)
- ~~architecture-extract-mugration-domain-logic.md~~ (resolved)
- ~~architecture-extract-prune-domain-logic.md~~ (resolved)
- ~~architecture-move-shared-domain-enums-from-commands.md~~ (resolved)
- [architecture-move-commands-to-cli-crate.md](architecture-move-commands-to-cli-crate.md)
- [architecture-add-treetime-python-crate.md](architecture-add-treetime-python-crate.md)
- ~~architecture-add-treetime-desktop-crate.md~~ (resolved: scaffolded with napi-rs async wrappers)
