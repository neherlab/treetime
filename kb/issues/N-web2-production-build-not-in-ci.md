# The production build of web2 does not run in CI

The `CI` workflow builds the WebAssembly back end of web2 in the dev mode and runs its smoke test, but no workflow runs the production build of web2. A change in the shared Rust crates, in the dist profile, in the toolchain, in the dev image, or in the web2 page can break the shipped web2 app without a failing check.

## Evidence

- The `CI` workflow runs one job per check group on every pull request and on pushes to `rust`, each through `dev/docker/run just check-group <group>` [.github/workflows/cli.yml#L64-L65](../../.github/workflows/cli.yml#L64-L65)
- The `clippy` group runs `lint-wasm`, and the `tests` group runs `test-wasm` [justfile#L107-L109](../../justfile#L107-L109). `test-wasm` builds the module with `build-wasm "dev"` and runs the Bun smoke test `packages/app-wasm/smoke.ts` [justfile#L269-L270](../../justfile#L269-L270)
- `build-web2` builds the module with `build-wasm "prod"` (the dist profile, as shipped) and then the Vite page with `bun run build:web2` [justfile#L592-L593](../../justfile#L592-L593). No check group and no workflow in `.github/workflows/` runs it, and `test-ts` builds only the web app (`build:web`)
- The module depends on version pairs that must change together: rustc and Emscripten (the Emscripten ABI has no semantic versioning, see the platform support page of `wasm32-unknown-emscripten` in the rustc book), Emscripten and the prebuilt archives of OpenBLAS (`binarylandia/build_openblas`), liblzma, and libzstd (`binarylandia/staticbourg`), each of which links only with the Emscripten version that built it, and the `wasm-bindgen` crate and the wasm-bindgen CLI that `emcc` runs. The dev-mode build in `test-wasm` covers these pairs, but not the flags of the dist profile

## Fix direction

Add `build-web2` to a check group, so the production build of the module and the web2 page run on every pull request.

## Validation

- A pull request that breaks the production build of the WebAssembly back end or of the web2 page fails a CI check
