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
- `app-napi`: napi addon (cdylib), AsyncTask wrappers calling app-api, `define_task!` macro
- `app-server`: axum HTTP API, `define_handler!` macro, spawn_blocking for CPU-bound work, JSON error responses, CORS

npm packages:

- `@neherlab/app-contracts`: TreeTimeBridge interface, arg types, CommandResult, ProgressEvent
- `@neherlab/app-ui`: shared React components, BridgeProvider/useBridge for dependency injection
- `@neherlab/app-desktop`: Electron main (IPC handlers loading napi addon) + preload (contextBridge) + Vite renderer
- `@neherlab/app-web`: Vite SPA, HTTP bridge via fetch()

## Bridge pattern

The React UI calls `bridge.ancestral(args)` and gets `Promise<CommandResult>`. The bridge implementation is injected via React context:

- Desktop: preload exposes IPC bridge via contextBridge -> ipcMain.handle -> napi addon AsyncTask -> Rust
- Web: fetch() -> Vite proxy -> app-server HTTP API -> tokio::spawn_blocking -> Rust

Both bridges return `CommandResult` (`{status: "ok"}` or `{status: "error", error: "..."}`) on all paths.

## Server auto-restart

Bacon watches Rust source files and restarts the server on changes. Config in `bacon.toml` with `on_change_strategy = "kill_then_restart"`. Uses the project's `.build/docker` target directory for cached builds.

## Desktop dev workflow

The desktop dev script (`scripts/dev.mjs`) starts a Vite dev server for the renderer (port 5174), then launches Electron with `VITE_DEV_SERVER_URL` pointing at Vite. In production, Electron loads built files from `dist/renderer/index.html`.

## Progress streaming (prepared, not yet wired)

`app-api` defines `ProgressSink` trait with `NoopProgress` and `StderrProgress` implementations. Planned:

- Desktop: napi ThreadsafeFunction callback -> IPC -> renderer
- Server: WebSocket or SSE -> browser
- CLI: StderrProgress (exists but not used by CLI yet)

## Open work

- Wire progress streaming implementations (napi ThreadsafeFunction, SSE/WebSocket)
- Extract I/O from algorithm library (run\_\* functions currently write to disk)
- Add structured result types (return data instead of writing files)
- Configure oxlint and prettier (packages installed, config files missing)
- Add vitest for frontend unit tests
- TypeScript arg type codegen or validation against Rust structs
