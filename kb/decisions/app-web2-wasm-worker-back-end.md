# Web2 runs the API router as WebAssembly in a web worker

**Type**: Application architecture without a v0 counterpart.

**Status**: Accepted. The back end builds for `wasm32-unknown-emscripten` and computes runs; the remaining work is in [Issues](#issues).

## Decision

The web2 app (`packages/app-web2`) renders the interface of `packages/app-ui` unchanged and reaches the back end through the same generated API client (`@neherlab/app-contracts/client`) as the web app and the desktop app. Only the `fetch` of the client differs:

- **Page**: `createPortFetch()` of `@neherlab/app-ui/port-fetch`, the fetch that the desktop renderer uses, sends each request through a message port. `WorkerConnection` (`packages/app-web2/src/backend/connection.ts`) compiles the WebAssembly module once, starts a dedicated worker, and hands the worker the compiled module and its end of a message channel
- **Worker**: the worker creates an instance of the module from the compiled module (the `instantiateWasm` hook of Emscripten, `instantiateBackend()` in `packages/app-web2/src/backend/instantiate.ts`), and a rejection of the hook fails the start. `ExchangeServer` (`packages/app-web2/src/backend/serve.ts`) parses each port message, calls `Backend.fetch()` of the instance, and pulls the response head and the body chunks of the exchange, which it sends back with the chunk buffers transferred. An abort message aborts the exchange
- **WebAssembly**: `Backend` (`packages/app-wasm/src/backend.rs`) holds the router of `api_router()`, the router of the public web server, so web2 has the routes and the validation of `treetime-server` and none of the routes of a local app. `Exchange` (`packages/app-wasm/src/exchange.rs`) drives one request through the router and gives its head and its body chunks on request, so a reader that stops pulling stops the body, and an abort ends a waiting read
- **Restart**: a worker that fails before it is ready stops the back end for good (`restarts: false`), because a restart would fail the same way. A worker that crashes later, or whose module aborts after a panic, is replaced by a new worker with the same compiled module (`restarts: true`); requests in flight fail and later requests wait for the new worker, as in the desktop app

## Target and platform

`app-wasm` is a binary crate for `wasm32-unknown-emscripten`. rustc links it through `emcc`, and `emcc` runs the wasm-bindgen CLI on the linked module (`-sWASM_BINDGEN`, experimental in wasm-bindgen 0.2.129). The result is one ES module, with the Emscripten runtime and the bindings, whose default export creates an instance, and the WebAssembly binary. The link settings are in `[target.wasm32-unknown-emscripten]` of `.cargo/config.toml`.

Emscripten provides the platform that the back end uses on the other targets:

- **C library**: OpenBLAS, LAPACK, and the compression libraries are compiled with `emcc` against the C library of Emscripten. The dev image installs a single-threaded OpenBLAS for the `WASM128_GENERIC` target, with LAPACK from the C translation that ships with OpenBLAS, which [binarylandia/build_openblas](https://github.com/binarylandia/build_openblas) builds with `emcc` for the Emscripten version of the image (`dev/docker/files/install-openblas-emscripten`), so the linear algebra of web2 goes through the same BLAS and LAPACK interface as the CLI. The functions of the C translation take no hidden string lengths and return `void`, as the `lapack-sys` declarations do; a call whose signature differs from its callee traps in WebAssembly. `xz` goes through `liblzma`, `zstd` through `zstd-sys`, both linked from the prebuilt WebAssembly archives of `binarylandia/staticbourg`
- **Files**: `std::fs` works on the in-memory file system of Emscripten, so the run store, the uploads, and the inputs and outputs of the commands run unchanged under `/treetime/runs`. The examples folder `/treetime/examples` is empty
- **Clocks, process id, entropy**: `SystemTime`, `Instant`, `process::id()`, and `getrandom` (through `getentropy`) work, so a run without `--seed` draws and logs its seed as on the other targets
- **Memory**: the module grows its memory up to the 4 GiB of `wasm32` and has a stack of 8 MiB

The module has no threads, no tokio runtime, and no unwinding. The back end adapts to each with one `cfg(target_family = "wasm")` alternative:

- **Threads**: a run computes in the calling thread instead of a thread of its own (`execute()` in `packages/app-backend/src/bridge/service.rs`). rayon runs its work on the calling thread, because `pthread_create` returns `ENOTSUP`, which Rust reports as `ErrorKind::Unsupported`
- **tokio runtime**: the blocking work of a request runs in place, where natively it runs on the blocking pool of tokio (`run_blocking()` and `read_body_blocking()` in `packages/app-server/src/api/blocking.rs`). A request body is read in full before the blocking work. A run file or a run archive is sent in chunks of at most 1 MiB that the exchange pulls one at a time (`FileChunks` and `RunZipChunks` in `packages/app-backend/src/runs/files.rs`, used by `packages/app-server/src/api/download.rs`), so a download never holds a whole output in the linear memory, where one allocation of 2 GiB or more is impossible in `wasm32`. The request timeout (`packages/app-server/src/api/timeout.rs`) and the keep-alive of the event streams need tokio timers and exist only in the native builds; both protect sockets, and the worker has none
- **Blocking HTTP client**: the examples download of the local apps (`packages/app-commands/src/examples_fetch.rs`) exists only in the native builds, because `reqwest` has no blocking client on WebAssembly targets
- **Unwinding**: wasm-bindgen does not support `panic = "unwind"` on this target, so a panic aborts the module, and the page restarts the worker

## Types at each boundary

- **API requests and responses**: the generated types and zod schemas of the OpenAPI document, as in the other apps
- **JavaScript interface of the module**: wasm-bindgen generates it from the Rust signatures, and `emcc --emit-tsd` writes its declarations (`just gen wasm-types` writes `packages/app-wasm/pkg/app_wasm.d.ts`). The worker takes the types of `Backend`, `Exchange`, and `ResponseHead` from the declared `MainModule`. The interface uses only primitives, byte arrays, and exported classes, so it needs no serde derives (the workspace serializes with deser) and no hand-written declaration. The one string type is the callback of `setPanicSink`, `(text: string) => void`, because wasm-bindgen declares a typed `js_sys::Function` as `Function`. The one hand-written type is the `instantiateWasm` option of the Emscripten runtime, which the declarations type as `unknown`
- **Message port messages** (`PortRequest`, `PortMessage`, `PortReply`): zod schemas declared once in `packages/app-ui/src/port.ts`. Every reply that the page receives and every message that the worker receives is parsed; a reply that breaks the protocol fails the open exchanges and later requests until a new port connects. These messages pass between two JavaScript contexts and carry the HTTP bytes as opaque data, so this is an exception to the rule that a value that reaches Rust is a Rust type with a derived schema: in web2 the worker turns them into the primitive arguments of the generated interface, and in the desktop app a type test (`packages/app-desktop/src/__tests__/test_port_types.test.ts`) keeps the types that napi-rs generates for the addon within these schemas
- **Worker control messages** (start, ready, failed, crashed): zod schemas in `packages/app-web2/src/backend/protocol.ts`; they never reach Rust

## Reason

- One implementation of each API route, check, and command for the CLI, the web app, the desktop app, and web2; the interface needs no second client
- The transport of the desktop app already sends API requests through a message port to an in-process router, so web2 reuses its fetch and its restart behavior instead of a second protocol
- An RPC interface per operation (the approach of treeknit and Nextclade, with comlink and tsify) would need a second client in `app-ui` for every operation, and tsify needs serde derives on the types
- Emscripten supplies the C library, the file system, the clocks, and the entropy that the back end uses on the other targets, so web2 keeps OpenBLAS and the C codecs, and the run store and the commands run without a file system abstraction. `wasm32-unknown-unknown` has none of them: it would need a C library written by the project for OpenBLAS, a file system abstraction in the I/O layer of `treetime-utils` and the run store, and browser back ends for clocks and entropy. WASI targets have no wasm-bindgen support, so the interface of the module would be hand-written

## Implementation

- `packages/app-wasm/src/backend.rs`, `exchange.rs`, `start.rs`: the WebAssembly back end; `packages/app-wasm/smoke.ts`: the smoke test of the built module that `just test-wasm` runs in Bun
- `packages/app-server/src/port.rs`: the requests and response headers of a message port exchange, shared with `packages/app-napi`
- `packages/app-server/src/api/blocking.rs`, `download.rs`, `timeout.rs`: the request work with and without a tokio runtime
- `packages/app-ui/src/port.ts`, `packages/app-ui/src/port-fetch.ts`: the message port messages and the fetch over a port
- `packages/app-web2/src/backend/`: the worker connection, the worker, and the exchange server
- `packages/app-ui/build/content-security-policy-vite.ts`: `'wasm-unsafe-eval'` for the page of web2, which compiles the module
- `dev/docker/files/install-emscripten`, `install-openblas-emscripten`: Emscripten and the prebuilt OpenBLAS in the dev image

## Issues

- [kb/issues/M-web2-run-blocks-the-back-end-worker.md](../issues/M-web2-run-blocks-the-back-end-worker.md)
- [kb/issues/M-web2-runs-are-lost-on-reload.md](../issues/M-web2-runs-are-lost-on-reload.md)
- [kb/issues/M-web2-offers-no-example-datasets.md](../issues/M-web2-offers-no-example-datasets.md)
- [kb/issues/M-web2-wasm-panic-loses-the-back-end-state.md](../issues/M-web2-wasm-panic-loses-the-back-end-state.md)
- [kb/issues/M-web2-supported-dataset-size-undocumented.md](../issues/M-web2-supported-dataset-size-undocumented.md)
- [kb/issues/M-web2-not-built-or-tested-in-ci.md](../issues/M-web2-not-built-or-tested-in-ci.md)
- [kb/issues/M-web2-not-verified-in-browsers.md](../issues/M-web2-not-verified-in-browsers.md)
- [kb/issues/M-web2-depends-on-experimental-wasm-bindgen-emscripten-support.md](../issues/M-web2-depends-on-experimental-wasm-bindgen-emscripten-support.md)
- [kb/issues/H-timetree-outputs-depend-on-the-rounding-of-the-platform.md](../issues/H-timetree-outputs-depend-on-the-rounding-of-the-platform.md)
