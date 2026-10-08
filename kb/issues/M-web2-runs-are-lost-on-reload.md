# The runs of web2 are lost when the page reloads

The run store of the WebAssembly back end lives in the in-memory file system of Emscripten under `/treetime/runs`. The run list, the run comparison, and the downloads work while the page is open; a reload of the page, or a restart of the worker after a panic ([M-web2-wasm-panic-loses-the-back-end-state.md](M-web2-wasm-panic-loses-the-back-end-state.md)), starts with an empty store.

## Evidence

- `browser_config()` in `packages/app-wasm/src/backend.rs` names `/treetime/runs`; `RunManager::open()` creates it in the file system of the instance
- No file system of Emscripten other than the in-memory one is mounted; the module links no IDBFS or WasmFS back end

> [!IMPORTANT]
> **Decision required.** Whether web2 keeps runs across page loads, and where:
>
> - no persistence: web2 holds the runs of one page session
> - IDBFS of Emscripten at `/treetime/runs`, synchronized to IndexedDB with `FS.syncfs()` after each change of the store
> - the WasmFS OPFS back end at `/treetime/runs`, which keeps the files in the origin private file system of the browser and which a job worker of [M-web2-run-blocks-the-back-end-worker.md](M-web2-run-blocks-the-back-end-worker.md) could share. The back end calls OPFS through a proxy worker thread, so it needs Emscripten pthreads, or JSPI without them (`system/lib/wasmfs/backends/opfs_backend.cpp` of Emscripten 6.0.10); the module has neither, and the JSPI support of wasm-bindgen on Emscripten needs Emscripten lifecycle hooks that no Emscripten release has (the "JSPI" section of the wasm-bindgen guide, `guide/src/reference/emscripten.md`)
>
> Browsers may evict both stores under storage pressure unless the page requests persistent storage.

## Validation

- After a reload of the page, the run list shows the runs of the previous session and their files download
