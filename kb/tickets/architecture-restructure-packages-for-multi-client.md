# Restructure existing packages for multi-client architecture

## Description

Reorganize the existing workspace so that the `treetime` library crate is a pure domain library consumable by any client (CLI, Python, desktop, WASM). The `treetime-cli` crate becomes one of several thin client crates.

This is an umbrella ticket coordinating prerequisite work items. Each item has its own ticket.

## Execution order

The work items have a dependency chain. Execute in this order:

### 1. Remove clap from library domain types

Ticket: `architecture-remove-clap-from-domain-types.md`

Remove `clap` from `treetime` crate's `[dependencies]`. Move CLI-specific derives to wrapper types in `treetime-cli` or `commands/*/args.rs`. Domain types use `strum` for string parsing.

### 2. Extract remaining domain logic from commands/

Three extraction tickets, independent of each other but all depend on step 1:

- `architecture-extract-coalescent-from-timetree.md` - move `commands/timetree/coalescent/` (1444 lines) to `src/coalescent/`
- `architecture-extract-timetree-inference-from-commands.md` - move inference (479 lines) and optimization (802 lines) to `src/timetree/`
- `architecture-extract-mugration-domain-logic.md` - move GTR refinement (381 lines) and discrete marginal (170 lines) to `src/mugration/`

### 3. Break dependency cycles in library

Two tickets, independent of each other, depend on step 2:

- `architecture-break-representation-gtr-cycle.md` - parameterize GTR inference over iterators
- `architecture-break-representation-ancestral-cycle.md` - move pure functions to `seq/`

### 4. Move commands/ shell to CLI crate

After steps 1-3, `commands/` contains only thin orchestration (`args.rs`, `run.rs`, `output/`). Move the entire `commands/` module to `treetime-cli`. The library crate no longer has a `commands/` module.

No existing ticket for this step - create when prerequisites are met.

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
  +-- representation/
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
  |   +-- shared/
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

- `treetime` crate compiles without `clap` dependency
- `treetime` crate has no `commands/` module
- `treetime-cli` builds and passes all existing tests
- `treetime-python` builds a wheel and runs a smoke test
- `treetime-desktop` compiles (frontend deferred)

## Related issues

Source: [H-core-multi-client-architecture-library-purity](../issues/H-core-multi-client-architecture-library-purity.md)

Coordinates:

- [architecture-remove-clap-from-domain-types.md](architecture-remove-clap-from-domain-types.md)
- [architecture-extract-coalescent-from-timetree.md](architecture-extract-coalescent-from-timetree.md)
- [architecture-extract-timetree-inference-from-commands.md](architecture-extract-timetree-inference-from-commands.md)
- [architecture-extract-mugration-domain-logic.md](architecture-extract-mugration-domain-logic.md)
- [architecture-break-representation-gtr-cycle.md](architecture-break-representation-gtr-cycle.md)
- [architecture-break-representation-ancestral-cycle.md](architecture-break-representation-ancestral-cycle.md)
- [architecture-add-treetime-python-crate.md](architecture-add-treetime-python-crate.md)
- [architecture-add-treetime-desktop-crate.md](architecture-add-treetime-desktop-crate.md)
