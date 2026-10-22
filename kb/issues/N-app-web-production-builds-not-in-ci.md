# The production build of web2 does not run in CI

The `CI` workflow lints both WebAssembly modules of web2 and runs the smoke test of the single-threaded module in the dev mode, but no workflow runs the production build of web2 or the smoke test of the threaded module. A change in the shared Rust crates, in the dist profile, in the toolchain, in the dev image, or in the web2 page can break the shipped web2 app without a failing check.

## Evidence

- The `CI` workflow runs one job per check group on every pull request and on pushes to `rust`, each through `dev/docker/run just check-group <group>` [.github/workflows/cli.yml#L64-L65](../../.github/workflows/cli.yml#L64-L65)
- The `clippy` group runs `lint-wasm` and `lint-wasm-mt`, and the `tests` group runs `test-wasm` [justfile#L118-L120](../../justfile#L118-L120). `test-wasm` builds the single-threaded module with `build-wasm "dev"` and runs the Bun smoke test `packages/app-wasm/smoke.ts` [justfile#L298-L299](../../justfile#L298-L299). `test-wasm-mt`, which runs the same smoke test under Node on the threaded module, is outside every gate [justfile#L301-L309](../../justfile#L301-L309)
- `build-web2` builds both modules with `build-wasm "prod"` and `build-wasm-mt "prod"` (the dist profile, as shipped) and then the Vite page with `bun run build:web2` [justfile#L641-L642](../../justfile#L641-L642). No check group and no workflow in `.github/workflows/` runs it, and `test-ts` builds only the web app (`build:web`)
- The module depends on version pairs that must change together: rustc and Emscripten (the Emscripten ABI has no semantic versioning, see the platform support page of `wasm32-unknown-emscripten` in the rustc book), Emscripten and the prebuilt archives of OpenBLAS (`binarylandia/build_openblas`), liblzma, and libzstd (`binarylandia/staticbourg`), each of which links only with the Emscripten version that built it, the `wasm-bindgen` crate and the wasm-bindgen CLI that `emcc` runs, and the nightly of the threaded build and the threaded OpenBLAS archive. The dev-mode build in `test-wasm` covers the pairs of the single-threaded module, and `lint-wasm-mt` compiles the threaded module without linking it, but neither covers the flags of the dist profile

## Fix direction

Add `build-web2` and `test-wasm-mt` to a check group, so the production build of both modules, the threaded smoke test, and the web2 page run on every pull request.

## Validation

- A pull request that breaks the production build of the WebAssembly back end or of the web2 page fails a CI check, and so does a pull request that breaks a run of the threaded module
