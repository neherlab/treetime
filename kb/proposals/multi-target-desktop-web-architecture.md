# Multi-target architecture: desktop and web apps sharing UI and backend

## Motivation

TreeTime needs a desktop app (replacing the legacy Flask web app) and a web app (for browser access). Both share the same React UI and call the same Rust algorithms. The architecture separates concerns into shared layers and target-specific shells.

## Architecture

```
treetime (algorithms)
    |
app-api (service layer: params, progress, commands)
    |
+---+----------+
|              |
desktop      server
(napi)       (axum)
```

TypeScript packages:

```
app-contracts (bridge interface, arg/result types)
    |
app-ui (shared React components, BridgeContext)
    |
+---+-----------+
|               |
desktop-app   web-app
(Electron)    (Vite SPA)
```

## Package inventory

Rust crates:

- `app-api`: ProgressSink trait, command wrappers accepting args + progress, re-exported arg types
- `app-napi`: napi addon (cdylib), AsyncTask wrappers calling app-api
- `app-server`: axum HTTP API, spawn_blocking for CPU-bound work, CORS

npm packages:

- `@neherlab/app-contracts`: TreeTimeBridge interface, arg types, result types, progress events
- `@neherlab/app-ui`: shared React components, BridgeProvider/useBridge for dependency injection
- `@neherlab/app-desktop`: Electron main + preload with contextBridge IPC
- `@neherlab/app-web`: Vite + React entry, HTTP bridge via fetch()

## Bridge pattern

The React UI calls `bridge.ancestral(args)` and gets `Promise<CommandResult>`. The bridge implementation is injected via React context:

- Desktop: preload exposes napi addon via contextBridge -> ipcMain -> napi AsyncTask -> Rust
- Web: fetch() -> app-server HTTP API -> tokio::spawn_blocking -> Rust

## Progress streaming (prepared, not yet wired)

`app-api` defines `ProgressSink` trait. Implementations:

- CLI: stderr output
- Desktop: napi ThreadsafeFunction callback -> IPC -> renderer
- Server: WebSocket or SSE -> browser

## Open work

- Wire napi addon IPC handlers in desktop main process
- Add `treetime serve` subcommand to CLI
- Extract I/O from algorithm library (run\_\* functions currently write to disk)
- Add structured result types (return data instead of writing files)
- Set up bun workspaces for npm packages
- Add oxlint + oxfmt for TypeScript linting/formatting
- Add Vite config for web-app dev server
- Add vitest for frontend unit tests
