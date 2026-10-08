# A panic in the web2 back end loses the state of the back end

The WebAssembly back end is built with `panic = "abort"`: a panic aborts the module, and the module cannot run again. The page restarts the worker with a new instance, so every run that the instance kept in its in-memory file system is gone.

## Evidence

- `-Cpanic=abort` in `[target.wasm32-unknown-emscripten]` of `.cargo/config.toml`: wasm-bindgen does not support unwinding across its boundary on Emscripten ([wasm-bindgen#5165](https://github.com/wasm-bindgen/wasm-bindgen/issues/5165)), although Emscripten itself can unwind with WebAssembly exceptions
- `CatchPanicLayer` in `finish_router()` of `packages/app-server/src/routes.rs` turns a panic into a 500 response with `catch_unwind`, which does not catch an abort
- `packages/app-wasm/src/start.rs` sends the panic text to the sink of `setPanicSink`; `packages/app-web2/src/backend/worker.ts` reports it as `crashed`, and `WorkerConnection` restarts the worker

> [!IMPORTANT]
> **Decision required.** Accept the abort and restart, and limit the loss with the persistence of [M-web2-runs-are-lost-on-reload.md](M-web2-runs-are-lost-on-reload.md), or build with `panic = "unwind"` once wasm-bindgen supports it on Emscripten.

## Validation

- After a forced panic in a route, the page reports the panic text, the worker restarts, and the runs of the persistent store are listed again
