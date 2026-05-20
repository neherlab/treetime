# Add treetime-desktop crate as thin library wrapper with Electron frontend

## Description

Create `packages/treetime-desktop/` as a workspace crate providing a minimal desktop application wrapping the Rust library. The Rust side exposes an IPC/command API. The frontend is a minimal Electron (or Tauri) shell.

The crate serves as structural backpressure: a second non-CLI consumer that depends on the `treetime` library crate, reinforcing that domain logic must stay in the library. Implementation of the full desktop UI is deferred - the initial crate establishes the architecture and a minimal proof-of-concept.

## Crate structure

```
packages/treetime-desktop/
  Cargo.toml          # workspace inheritance, tauri or napi dependencies
  src/
    lib.rs            # command/IPC handlers calling treetime library
    commands.rs       # desktop command definitions (run_ancestral, run_clock, etc.)
  frontend/           # minimal Electron/Tauri frontend (deferred)
    package.json
    src/
      index.html
      main.js
```

## Dependencies

- `treetime` (workspace) - core library
- `treetime-io` (workspace) - file format I/O
- `treetime-utils` (workspace) - error handling
- `serde` + `serde_json` (workspace) - IPC serialization
- `tauri` or `napi-rs` - desktop integration (decision deferred)

## Scope

Initial implementation:

- Cargo.toml with workspace inheritance
- `src/lib.rs` with a single command handler calling `treetime::commands::ancestral`
- No frontend UI (placeholder HTML only)

Full desktop UI implementation is out of scope for this ticket.

## Prerequisites

- `commands/` module moved to `treetime-cli` (ticket: `architecture-move-commands-to-cli-crate.md`)

## Related issues

Source: [H-core-multi-client-architecture-library-purity](../issues/H-core-multi-client-architecture-library-purity.md)
