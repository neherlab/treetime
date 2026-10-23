# The production build of the web app and the threaded smoke test do not run in CI

The `Web app` workflow builds web2 for production on every push, but no check runs the production build of the web app (`app-web`, with the API server) or the smoke test of the threaded WebAssembly module. A change in a Vite config or build plugin, in a page, or in the threaded module can break the shipped web app or a threaded run without a failing check.

## Evidence

- The `CI` workflow runs one job per check group on every pull request and on pushes to `rust`, each through `dev/docker/run just check-group <group>` [.github/workflows/cli.yml#L65](../../.github/workflows/cli.yml#L65)
- The `clippy` group runs `lint-wasm` and `lint-wasm-mt`, and the `tests` group runs `test-wasm` [justfile#L124-L126](../../justfile#L124-L126). `test-wasm` builds the single-threaded module with `build-wasm "dev"` and runs the Bun smoke test `packages/app-wasm/smoke.ts` [justfile#L301-L302](../../justfile#L301-L302). `test-wasm-mt`, which runs the same smoke test under Node on the threaded module, is outside every gate [justfile#L307](../../justfile#L307)
- The `Web app` workflow runs `just build-web2` on every push to any branch and when a pull request opens or reopens, in the image of the WebAssembly target [.github/workflows/web.yml#L65](../../.github/workflows/web.yml#L65). It builds both modules in the dist profile, as shipped, and the Vite page, so a commit that breaks the production build of web2 gets a failing `Web app` check
- `build-web` builds the web app with `bun run build:web` [justfile#L649-L650](../../justfile#L649-L650). No check group and no workflow runs it. `test-ts` runs `bun test --parallel` and builds no app
- Vite loads each config in Node and leaves workspace imports of the config to Node. A config that Node cannot load breaks every build, dev server, and image of that app, while `just check` and `just test-ts` pass: the example-data plugin once imported the extensionless entry of `@neherlab/app-contracts`, and only a manual `just build-web2` showed it
- The threaded module depends on version pairs that must change together: the nightly of the threaded build and its `-Zbuild-std`, Emscripten and the threaded OpenBLAS archive (`binarylandia/build_openblas`), which links only with the Emscripten version that built it, and the pthread pool of the module. The production build of the `Web app` workflow links the threaded module but never runs it

## Fix direction

Add `build-web` and `test-wasm-mt` to a check group, so the production build of the web app and a run of the threaded module happen on every pull request.

## Validation

- A pull request that breaks the production build of the web app, or a run of the threaded module, fails a CI check
