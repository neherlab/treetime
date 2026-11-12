---
title: Web2 app
description: The web app on the WebAssembly back end in a web worker.
---

The web2 app (`packages/app-web2`) is the web app without a server: the interface of `app-ui` sends its API requests to the WebAssembly build of the API router (`packages/app-wasm`), which runs in a dedicated web worker of the page. [kb/decisions/app-web2-wasm-worker-back-end.md](https://github.com/neherlab/treetime/blob/rust/kb/decisions/app-web2-wasm-worker-back-end.md) describes the design.

| Task                                                                       | Command                                   |
| -------------------------------------------------------------------------- | ----------------------------------------- |
| Build the single-threaded WebAssembly module into `packages/app-wasm/pkg/st/` | `just build-wasm [dev\|release\|prod]`    |
| Build the threaded WebAssembly module into `packages/app-wasm/pkg/mt/`     | `just build-wasm-mt [dev\|release\|prod]` |
| Clippy of the single-threaded and of the threaded module                   | `just lint-wasm`, `just lint-wasm-mt`     |
| Smoke test of the single-threaded module, in Bun                           | `just test-wasm`                          |
| Smoke test of the threaded module with 2 and 8 threads, in Node            | `just test-wasm-mt`                       |
| Start the web2 dev server                                                  | `just up-web2`                            |
| Check that the web2 dev server answers and runs the checkout HEAD           | `just health`                             |
| Build web2 for production and serve it                                     | `just serve-web2`                         |
| Production build of web2 into `packages/app-web2/dist/`                    | `just build-web2`                         |

`just build-wasm` builds `app-wasm` for `wasm32-unknown-emscripten` with cargo. rustc links the module through `emcc` of Emscripten, and `emcc` runs the `wasm-bindgen` CLI on the linked module, so the build writes one JavaScript module with the Emscripten runtime and the bindings, the WebAssembly binary, and the TypeScript declarations; the recipe copies them into the package `@neherlab/app-wasm` (`packages/app-wasm/pkg/`). Emscripten provides the C library of OpenBLAS and of the compression libraries, an in-memory file system for the runs, clocks, and entropy. The dev image installs Emscripten, a single-threaded and a threaded OpenBLAS that [binarylandia/build_openblas](https://github.com/binarylandia/build_openblas) builds with `emcc` for this Emscripten version, and the WebAssembly liblzma and libzstd of staticbourg, all from `.config/mise.image.toml` (the `http:emscripten` and `*-wasm32-unknown-emscripten` entries). Emscripten gives no ABI guarantee between versions, so the image build fails when an archive comes from another Emscripten version, and an Emscripten upgrade needs new `build_openblas` and staticbourg releases first. The `prod` mode uses the `dist` profile, which `emcc` optimizes at link time. `just gen wasm-types` writes the declarations `packages/app-wasm/pkg/app_wasm.d.ts`, the only committed file of the package, so the TypeScript checks run without a WebAssembly build. The version of the `wasm-bindgen` CLI in `.config/mise.toml` and the version of the `wasm-bindgen` crate in `Cargo.toml` must be equal, because the CLI accepts only modules of its own schema version. mise installs the CLI from [binarylandia/appville](https://github.com/binarylandia/appville) (task `tasks/wasm-bindgen`): upstream wasm-bindgen with a patch that keeps the start function of the linker, which the threaded module needs to set up its shared memory. A new crate version needs a new appville release first.

The common link settings of both modules are in `[target.wasm32-unknown-emscripten]` of `.cargo/config.toml`, among them `-sSTACK_OVERFLOW_CHECK=2` and `-Wl,--stack-first`, which stop the module with "stack overflow" instead of letting an overflow overwrite data. Each cargo call for the target adds the settings of its build with `--config`: `dev/cargo/wasm-st.toml` for the single-threaded module and `dev/cargo/wasm-mt.toml` for the threaded one. `just build-wasm-mt` builds on the nightly of `dev/docker/files/wasm-threads-toolchain`, which rebuilds the standard library with atomics (`-Zbuild-std`), in its own target folder, and links with `-pthread` and the threaded OpenBLAS, liblzma, and libzstd. The module starts a pool of `jobs + 2` pthread workers with the module, because rayon waits for every thread of its pool and a new pthread worker starts only after its creator returns to the event loop. `just build-wasm-mt` fails when the threaded module declares other exports than `pkg/app_wasm.d.ts`. The package exports the single-threaded module as `@neherlab/app-wasm` and the threaded one as `@neherlab/app-wasm/threads`.

`just lint-wasm` and `just lint-wasm-mt` run clippy over the binary of `app-wasm` for the target, the second over the code that compiles only with atomics; `just check` runs both. The tests of `app-wasm` run natively, because a test harness needs unwinding and the module aborts on a panic. `just test-wasm` builds the single-threaded module and runs `packages/app-wasm/smoke.ts` in Bun: it creates a back-end instance and a job instance of the module and drives the job protocol without workers, with an `ancestral` run of `data/flu/h3n2/20` and its files and archive, a request answered while a job computes, a cancelled run, and an input that fails its checksum. `just test-wasm-mt` runs the same smoke test on the threaded module with 2 and with 8 threads under Node, because Bun does not start the pthread workers of Emscripten, and Node stays alive during a threaded run only through a timer that the test holds; it is slow and outside every gate. `just up-web2` builds both modules once; run `just build-wasm` and `just build-wasm-mt` again after a Rust change. `just up-web2` and `just serve-web2` run in the foreground until Ctrl-C. `just health` and `just status` report whether the dev server runs and which commit it was started from.

The page compiles the single-threaded module once and starts the back-end worker with it (`WorkerConnection`, `packages/app-web2/src/backend/connection.ts`). The worker answers each request through an exchange of the module and sends the response head and body chunks back through a message port (`ExchangeServer`, `packages/app-web2/src/backend/serve.ts`). The page uses the fetch of the desktop renderer for this port (`createPortFetch()` of `@neherlab/app-ui/port-fetch`). An upload body passes to the worker as a `Blob`, which a structured clone passes by reference, and the worker feeds it to the exchange in slices (`bodyChunks()` in `serve.ts`), so the module never holds a whole upload in one buffer; before the transfer, the page compares the file size with the upload allowance that remains for the run. A back-end worker that fails to start stops the back end; a back-end worker that crashes restarts with the same module, and the requests in flight and every run of the page are lost. The content security policy of the web2 page allows `'wasm-unsafe-eval'`, which WebAssembly compilation needs, and nothing else that evaluates code.

Each run computes in a job worker of its own, which the back-end worker starts for the run (`JobHost`, `packages/app-web2/src/backend/jobs.ts`, and `packages/app-web2/src/job/worker.ts`), so several runs can compute at the same time. The back-end worker sends the job its input files, writes each progress event of the job into the run's event log, so the page shows the progress live, and receives the output files. Stop terminates the job worker at once. A panic, an out-of-memory abort, or a stack overflow in a job fails only that run, with the stderr lines and the abort reason of the module as the error. The job worker runs the threaded module, 8 threads by default; the URL parameter `jobs` sets the number of threads of each run (`?jobs=4`), and `?jobs=1` selects the single-threaded module. Threads need a cross-origin isolated page in a secure context (HTTPS or `localhost`): the dev server and the preview server of `just serve-web2` send `Cross-Origin-Opener-Policy: same-origin` and `Cross-Origin-Embedder-Policy: require-corp` (`packages/app-web2/vite.config.ts`). Without isolation, without `SharedArrayBuffer`, or when the browser cannot create a shared memory of 4 GiB, and after a threaded job fails to start, the jobs run the single-threaded module, and the console says why.

The back end serves the public routes of `api_router()` (`packages/app-api`), the router that the desktop app also serves, without the routes of `local_api_routers()` that only the desktop app has (settings, workspace, examples download). The runs, the uploads, and the outputs live in the in-memory file system of the back-end worker. [kb/decisions/app-web2-wasm-worker-back-end.md](https://github.com/neherlab/treetime/blob/rust/kb/decisions/app-web2-wasm-worker-back-end.md) lists the open issues.

## Run model

Runs live in the tab: nothing stores a run in browser storage, and a reload or a closed tab ends every run of the tab. A run leaves the tab as a run package, which the user downloads and imports again in any tab or in the desktop app. Desktop keeps its runs in its runs folder with the same back-end code, so the two run histories differ only in their lifetime. [kb/decisions/app-web2-runs-per-tab.md](https://github.com/neherlab/treetime/blob/rust/kb/decisions/app-web2-runs-per-tab.md) records the decision and the options weighed.

The back end keeps one layout on both apps (`RunStore`, `packages/app-backend/src/runs/store.rs`), rooted at `/treetime/runs` in the in-memory file system of web2 and at the runs folder on desktop:

- `<id>/run.json`: the `RunRecord`: the configuration, the inputs with their size and SHA-256, the output files, the status, the warnings, the seed, `source_run_id`, the derivation group id, and the size of the run. Every path in the record is relative to the run folder; desktop inputs from outside the runs folder stay absolute
- `<id>/events.jsonl`: the progress events of the run, which the Log view reads
- `<id>/out/`: the output files
- `.inputs/<sha256>/<name>`: the uploaded inputs, stored once by content (`InputStore`) and shared by every run that reads them, so a variant keeps working after the delete of its source. `DELETE /api/runs/{id}` removes a run and the inputs that no other run reads, and emits one `AppChange::RunRemoved`
- `.staging/`: the import of a run package, which the back end moves to `<id>/` only when it is complete

The page keeps the delete of a run hidden for the Undo time of its notice and sends one delete request per run when the notice closes (`useRunDeletionStore`, `packages/app-ui/src/store/runDeletion.ts`). The form deletes the upload run of a draft that the user discards. While the tab holds runs, the page asks before it unloads, and after a crash of the back-end worker it shows a banner that every run of the tab is lost.

## Run package

A run package is a zip file with one folder `treetime-<id>/` (`export_run()`, `packages/app-backend/src/runs/package.rs`), written as a stream by `GET /api/runs/{id}/package`:

| Path                  | Content                                                                                      |
| --------------------- | -------------------------------------------------------------------------------------------- |
| `run.json`            | The run record, with every path relative to the package folder                               |
| `events.jsonl`        | The progress events of the run                                                               |
| `inputs/<name>`       | Each input file; `inputs/<n>/<name>` when two inputs share a name                            |
| `out/`                | The output files                                                                             |
| `run.yaml`            | The configuration of the run, with the paths of the package and the seed of the run          |
| `command.sh`          | The `treetime` command line that runs `run.yaml` in the package folder and writes into `out/` |

`POST /api/runs/import` (`import_run()`, `packages/app-backend/src/runs/package_import.rs`) takes the zip as the request body. It checks that `run.json` parses as a record of an ended run, that every path is relative and inside the package, that each input has the size and SHA-256 of the record, and that every output file of the record exists; it gives the run a new id when the id exists, and publishes it through `.staging/` with one `RunCreated` event. The desktop app imports the packages of web2, and web2 those of the desktop app. The same `run.yaml` and `command.sh` reproduce the run with the CLI. The page saves a package and the output files through `saveRunOutput()`, which streams to disk through `showSaveFilePicker()` where the browser has it and builds one Blob in the other browsers.

## Ports

Each checkout takes its ports from its port slot, a group of 10 ports in the block 14000-14999, so worktrees run side by side. The main checkout holds slot 0, whose ports start at 14000. A linked worktree claims the lowest slot that no other checkout of the repository holds and whose ports have no listener, and keeps it until the worktree is removed; the claims are files in the git common directory (`dev/lib/ports.sh`). `just status` prints the slot and its ports. `TREETIME_WEB2_PORT`, `TREETIME_WEB2_SERVE_PORT`, and the other `TREETIME_*_PORT` variables override single ports.

A route launcher on the host can publish the web2 dev server through a local HTTPS proxy, such as [portless](https://github.com/vercel-labs/portless), under a hostname such as `https://treetime.localhost`. The launcher starts `just up-web2` as its child process, assigns the hostnames, and removes the route when the command ends. The checkout and the launcher exchange two files in `tmp/portless/`:

- `web.hostnames.json`: `{"hostnames": [...], "urls": [...]}`, written by the launcher before it starts the command. `just status` prints the URLs. Without this file, the dev server has no public hostname
- `web.json`: `{"port": <n>}`, written by `just up-web2` after the dev server answers, through a temporary file and a rename, so the launcher never reads a partial file. `just up-web2` writes it also when no launcher runs

The checkout never reads or writes the state of the proxy itself. The Vite server accepts every `*.localhost` hostname, so each assigned name works without configuration. The HTTPS proxy also gives the page the secure context that the threaded module needs on a hostname other than `localhost`.

## Browser storage

The browser keeps three keys (`browserPreferencesStorage()`, `packages/app-ui/src/preferences/storage.ts`):

- `treetime-preferences` in `localStorage`: the width of the sidebar (`UiSettings.sidebar_width`), shared by every tab
- `treetime-draft` in `sessionStorage`: the unfinished form (`UiDraft`), so each tab keeps its own form through a reload. A draft that does not parse is dropped. An example input keeps its example path (`UiDraftOrigin` `example`), and the page fetches it again after a reload; an uploaded file belonged to the lost back end, and a notice asks the user to choose it again
- `theme` in `localStorage`: the theme of `next-themes`, which applies a change in one tab to the other tabs at once

The runs stay in the memory of the back-end worker and are not stored in the browser.
