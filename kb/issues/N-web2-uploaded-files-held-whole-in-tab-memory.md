# web2 holds each uploaded file whole in the memory of the tab

Uploads reach the WebAssembly back end in slices and the input store writes them in chunks, but the in-memory file system of Emscripten keeps each stored file whole in JavaScript memory of the back-end worker, and the job worker of a run receives its inputs into its own in-memory file system. An upload near the 2 GiB upload limit therefore needs about that much memory in the back-end worker, and again in the job worker while the run computes, next to the outputs of every run of the tab. Browsers limit the memory of a tab, so such an upload can end the tab with every run.

## Evidence

- Uploads stream from the page to the input store (`bodyChunks()` in `packages/app-web2/src/backend/serve.ts`, `upload_body()` in `packages/app-api/src/api/blocking.rs`, `InputStore::begin()` in `packages/app-backend/src/runs/input_store.rs`)
- The runs root of web2 is `/treetime/runs` in the MEMFS of Emscripten (`RUNS_DIR` in `packages/app-wasm/src/backend.rs`), which holds each file as one array
- The job host transfers the input files of a run into its job worker (`JobHost`, `packages/app-web2/src/backend/jobs.ts`)
- [N-web2-supported-dataset-size-undocumented.md](N-web2-supported-dataset-size-undocumented.md) measured the Node process of a module at 1.8 to 2.7 GiB for the large examples, because MEMFS keeps the files outside the linear memory

> [!IMPORTANT]
> **Investigation required.** Measure the memory of the tab in Chrome and Firefox after an upload of 1 GiB and 2 GiB and during a run that reads it, and find the largest upload that the browsers keep. Options to weigh afterwards: a lower upload limit in web2, or a file system that keeps uploads outside the memory of the tab (WasmFS with the origin private file system, which [kb/decisions/app-web2-runs-per-tab.md](../decisions/app-web2-runs-per-tab.md) did not take for runs).

## Validation

- The upload limit of web2 matches the measured capacity of the browsers, or large uploads no longer stay in tab memory
