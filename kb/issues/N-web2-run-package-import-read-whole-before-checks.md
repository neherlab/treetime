# web2 reads a run package import whole into memory before it checks it

The page sends a run package to `POST /api/runs/import` in slices, but the WebAssembly back end collects the whole body into one buffer before the import reads it. A run package can be larger than 1 GB (ancestral on `mpox/clade-ii/2000` with default outputs is about 1.18 GB), so the import needs the whole package in the linear memory of the back-end worker, next to the runs of the tab, before the checks of the package start. A package near 2 GiB cannot be allocated in `wasm32`, and the abort loses every run of the tab ([M-web2-wasm-panic-loses-the-back-end-state.md](M-web2-wasm-panic-loses-the-back-end-state.md)). The native builds stream the body.

## Evidence

- `read_body_blocking()` on wasm: `to_bytes(body, usize::MAX)` [packages/app-api/src/api/blocking.rs#L41-L49](../../packages/app-api/src/api/blocking.rs#L41-L49); the native variant streams through `SyncIoBridge` (L30-L39)
- The import route calls it: `read_body_blocking(body, move |reader| runs.import(reader))` [packages/app-api/src/routes.rs#L663](../../packages/app-api/src/routes.rs#L663)
- Uploads stream on wasm through `upload_body()` in the same file, so only the import keeps the whole-body read

## Fix direction

Feed the body of the import to the zip reader in chunks on wasm, as `upload_body()` does for uploads: write the body to a temporary file of the runs root through the chunked writer, then open that file for the checks, or check a size limit before the read.

## Validation

- A Rust test imports a package through the wasm body path with a body of several chunks; the wasm smoke test imports a package of several chunks
