# A crash of the web2 back-end worker loses every run

The WebAssembly back end is built with `panic = "abort"`: a panic aborts the module, and the module cannot run again. A run computes in a job worker of its own, so a panic, an out-of-memory abort, or a stack overflow in a job ends only that run ([kb/decisions/app-web2-wasm-worker-back-end.md](../decisions/app-web2-wasm-worker-back-end.md)). A crash of the back-end worker itself still loses every run: the page terminates the worker, which also ends the job workers that it started, and restarts the back end with a new instance, whose in-memory file system holds no runs, no uploads, and no outputs.

## Evidence

- `-Cpanic=abort` in `[target.wasm32-unknown-emscripten]` of `.cargo/config.toml`: wasm-bindgen does not support unwinding across its boundary on Emscripten ([wasm-bindgen#5165](https://github.com/wasm-bindgen/wasm-bindgen/issues/5165)), although Emscripten itself can unwind with WebAssembly exceptions
- `CatchPanicLayer` in `finish_router()` of `packages/app-server/src/routes.rs` turns a panic into a 500 response with `catch_unwind`, which does not catch an abort
- The back-end worker reports an abort or a trap of its module as `crashed` with the error lines of `ErrorLines` (`packages/app-web2/src/backend/instantiate.ts`), and `WorkerConnection` (`packages/app-web2/src/backend/connection.ts`) terminates the worker and starts a new one; the job host of the old worker ends with it
- Causes that remain in the back-end worker: a panic in a route or in the preparation of a run, an allocation that fails while the worker holds the uploads and the outputs of every run of the page in its linear memory and its file system, and a trap in a call of the job host into the module

> [!IMPORTANT]
> **Decision required.** Accept the abort and restart, and limit the loss with a store that outlives the back-end worker (IndexedDB through IDBFS, or the origin private file system through WasmFS), or build with `panic = "unwind"` once wasm-bindgen supports it on Emscripten.

## Validation

- After a forced panic in a route, the page reports the panic text, the worker restarts, and the runs of the page are listed again, or the decision file states that a crash of the back-end worker loses them
