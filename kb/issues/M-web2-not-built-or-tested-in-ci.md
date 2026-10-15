# Web2 is not built or tested in CI

No workflow builds the WebAssembly back end of web2 or runs its smoke test, so a change in the shared Rust crates, in the toolchain, or in the dev image can break web2 without a failing check.

## Evidence

- `.github/workflows/` (`cli-build.yml`, `cli.yml`, `desktop-build.yml`, `nightly.yml`, `web-deploy.yml`) contains no `wasm32-unknown-emscripten` build and no `just test-wasm`
- `just check` lints the module (`just lint-wasm`) but builds no linked module; `just test-wasm` builds and runs it, and only a developer runs it
- The module depends on version pairs that must change together: rustc and Emscripten (the Emscripten ABI has no semantic versioning, see the platform support page of `wasm32-unknown-emscripten` in the rustc book), Emscripten and the prebuilt OpenBLAS archive of `binarylandia/build_openblas` (each archive links only with the Emscripten version that built it), and the `wasm-bindgen` crate and the wasm-bindgen CLI that `emcc` runs

## Fix direction

A CI job that builds the dev image layer with Emscripten and OpenBLAS, runs `just test-wasm`, and builds web2 with `just build-web2`.

## Validation

- A pull request that breaks the build or the smoke test of the WebAssembly back end fails a CI check
