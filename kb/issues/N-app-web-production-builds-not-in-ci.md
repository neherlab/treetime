# The production builds of the web apps do not run in CI

The `CI` workflow lints both WebAssembly modules of web2 and runs the smoke test of the single-threaded module in the dev mode, but no check runs the production build of the web app or of web2, or the smoke test of the threaded module. A change in the shared Rust crates, in the dist profile, in the toolchain, in the dev image, in a Vite config or build plugin, or in a page can break a shipped web app without a failing check.

## Evidence

- The `CI` workflow runs one job per check group on every pull request and on pushes to `rust`, each through `dev/docker/run just check-group <group>` [.github/workflows/cli.yml#L64-L65](../../.github/workflows/cli.yml#L64-L65)
- The `clippy` group runs `lint-wasm` and `lint-wasm-mt`, and the `tests` group runs `test-wasm` [justfile#L118-L120](../../justfile#L118-L120). `test-wasm` builds the single-threaded module with `build-wasm "dev"` and runs the Bun smoke test `packages/app-wasm/smoke.ts` [justfile#L298-L299](../../justfile#L298-L299). `test-wasm-mt`, which runs the same smoke test under Node on the threaded module, is outside every gate [justfile#L301-L309](../../justfile#L301-L309)
- `build-web2` builds both modules with `build-wasm "prod"` and `build-wasm-mt "prod"` (the dist profile, as shipped) and then the Vite page with `bun run build:web2` [justfile#L651-L652](../../justfile#L651-L652). No check group and no workflow in `.github/workflows/` runs it
- `build-web` builds the web app with `bun run build:web` [justfile#L646-L647](../../justfile#L646-L647). No check group runs it; only the deploy workflow builds the web image from a given commit [.github/workflows/web-deploy.yml#L62](../../.github/workflows/web-deploy.yml#L62). `test-ts` runs `bun test --parallel` and builds no app
- Vite loads each config in Node and leaves workspace imports of the config to Node. A config that Node cannot load breaks every build, dev server, and image of that app, while `just check` and `just test-ts` pass: the example-data plugin once imported the extensionless entry of `@neherlab/app-contracts`, and only a manual `just build-web2` showed it
- The module depends on version pairs that must change together: rustc and Emscripten (the Emscripten ABI has no semantic versioning, see the platform support page of `wasm32-unknown-emscripten` in the rustc book), Emscripten and the prebuilt archives of OpenBLAS (`binarylandia/build_openblas`), liblzma, and libzstd (`binarylandia/staticbourg`), each of which links only with the Emscripten version that built it, the `wasm-bindgen` crate and the wasm-bindgen CLI that `emcc` runs, and the nightly of the threaded build and the threaded OpenBLAS archive. The dev-mode build in `test-wasm` covers the pairs of the single-threaded module, and `lint-wasm-mt` compiles the threaded module without linking it, but neither covers the flags of the dist profile

## Fix direction

Add `build-web`, `build-web2`, and `test-wasm-mt` to a check group, so the production builds of both web apps and both modules, and the threaded smoke test, run on every pull request.

## Validation

- A pull request that breaks the production build of the WebAssembly back end, of the web app, or of web2 fails a CI check, and so does a pull request that breaks a run of the threaded module
