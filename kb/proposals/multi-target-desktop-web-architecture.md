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

`treetime` defines `ProgressSink` trait with `NoopProgress`. `app-api` has `StderrProgress`. Planned:

- Desktop: napi ThreadsafeFunction callback -> IPC -> renderer
- Server: SSE -> browser
- CLI: StderrProgress (exists but not used by CLI yet)

## Completed work

- `ProgressSink` trait moved from app-api to treetime core
- clap feature-gated in treetime (optional dependency behind `clap` feature flag)
- schemars dependency added, JSON Schema generation via `treetime schema write` CLI subcommand
- build.rs auto-generates JSON Schema into app-contracts on treetime source changes
- Result types defined for all 7 commands (run\_\* functions return structured data)
- Per-command TypeScript result types in app-contracts (AncestralResult, ClockResult, etc.)
- Per-command React hooks in app-ui (useAncestral, useClock, etc.)
- Desktop: fixed missing version IPC handler and QueryProvider
- Server: error response format changed to {code, message}
- Convert binary removed from app-cli
- Homoplasy stub returns HomoplasyResult type

## Open work

- Wire progress streaming implementations (napi ThreadsafeFunction, SSE)
- Extract I/O from run*\* into execute*\* functions (computation/I/O separation per command)
- Define computation input structs (parsed data, no file paths)
- Wire structured results through app-api, app-napi, app-server (currently discarded)
- Wire CLI through app-api (currently bypasses it)
- Set up TS codegen Turbo pipeline (json-schema-to-typescript + ts-to-zod)
- Add Serialize/JsonSchema to complex types (GTR, graph payloads) for full schema generation
- Desktop progress IPC, file dialog
- Web SSE progress, file upload
- Zustand state store in app-ui
- Add vitest for frontend unit tests
