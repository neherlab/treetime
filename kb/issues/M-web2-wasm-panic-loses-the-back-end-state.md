# A crash of the web2 back-end worker loses every run

The WebAssembly back end is built with `panic = "abort"`: a panic aborts the module, and the module cannot run again. A run computes in a job worker of its own, so a panic, an out-of-memory abort, or a stack overflow in a job ends only that run ([kb/decisions/app-web2-wasm-worker-back-end.md](../decisions/app-web2-wasm-worker-back-end.md)). A crash of the back-end worker itself still loses every run: the page terminates the worker, which also ends the job workers that it started, and restarts the back end with a new instance, whose in-memory file system holds no runs, no uploads, and no outputs. The page shows a persistent banner that every run of the tab is lost and that run packages keep runs; the runs of the tab are the expected loss of the run model of [kb/decisions/app-web2-runs-per-tab.md](../decisions/app-web2-runs-per-tab.md), in which the run package is the durable copy. What remains is the cause: a crash of the back-end worker still ends runs that the user has not downloaded yet, and nothing recovers them. After a failed start, the banner offers only "Reload the page", with no restart in place.

## Evidence

- `-Cpanic=abort` in `[target.wasm32-unknown-emscripten]` of `.cargo/config.toml`: wasm-bindgen does not support unwinding across its boundary on Emscripten ([wasm-bindgen#5165](https://github.com/wasm-bindgen/wasm-bindgen/issues/5165)), although Emscripten itself can unwind with WebAssembly exceptions
- `CatchPanicLayer` in `finish_router()` of `packages/app-api/src/routes.rs` turns a panic into a 500 response with `catch_unwind`, which does not catch an abort
- The back-end worker reports an abort or a trap of its module as `crashed` with the error lines of `ErrorLines` (`packages/app-web2/src/backend/instantiate.ts`), and `WorkerConnection` (`packages/app-web2/src/backend/connection.ts`) terminates the worker and starts a new one; the job host of the old worker ends with it
- Causes that remain in the back-end worker: a panic in a route or in the preparation of a run, an allocation that fails while the worker holds the uploads and the outputs of every run of the page in its in-memory file system, and a trap in a call of the job host into the module. Uploads reach the module in slices and no longer as one buffer, and `DELETE /api/runs/{id}` frees the files of a run, which lowers the memory that the runs of a tab hold

## Fix direction

Build with `panic = "unwind"` once wasm-bindgen supports it on Emscripten, so that `CatchPanicLayer` turns a panic in a route into an error response and the worker keeps its runs. A store that outlives the back-end worker is the option that [kb/decisions/app-web2-runs-per-tab.md](../decisions/app-web2-runs-per-tab.md) did not take.

## Validation

- After a forced panic in a route, the request fails with the panic text and the runs of the page stay listed
