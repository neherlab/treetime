# A run of web2 blocks the back-end worker until it ends

The WebAssembly back end has no threads, so `AppService::start()` computes a run in the calling thread (`execute()` in `packages/app-commands/src/bridge/service.rs`). The start request of a run returns when the run ends, and the back-end worker answers no other request in the meantime.

## Evidence

- `packages/app-commands/src/bridge/service.rs`: the `cfg(target_family = "wasm")` variant of `execute()` calls `finish()` directly; the native variant starts a thread per run
- Emscripten without pthreads returns `ENOTSUP` from `pthread_create`, so `std::thread::spawn` fails with `ErrorKind::Unsupported`; rayon then runs its work on the calling thread
- The page follows a run through its event stream (`GET /api/runs/{id}/events`). During the run the worker sends no events; afterwards the stream replays the events of the run from its event log
- The page itself stays responsive, because the run computes in the dedicated worker

> [!IMPORTANT]
> **Decision required.** How web2 computes runs beside the API:
>
> - one job worker per run, a second instance of the same module started by the page on request of the back-end worker. Each instance has its own in-memory file system, so the job worker needs the inputs of the run and returns its outputs: either copied as transferred buffers, or through one OPFS folder that both instances mount with the WasmFS OPFS back end of Emscripten, which needs pthreads or JSPI ([M-web2-runs-are-lost-on-reload.md](M-web2-runs-are-lost-on-reload.md))
> - WebAssembly threads: Emscripten pthreads need a Rust standard library built with `+atomics` (nightly Rust with `-Z build-std`) and `SharedArrayBuffer`, which needs cross-origin isolation (COOP and COEP headers) that static hosts such as GitHub Pages cannot send
> - keep the runs in the back-end worker, and show the page a waiting state while the start request is open

## Validation

- In web2, the page receives the progress events of a running run and can open another run while it computes
