# The production build of the web app and the threaded smoke test do not run in CI

The `Web app` workflow builds web2 for production on every push, but no check runs the production build of the web app (`app-web`, with the API server) or the smoke test of the threaded WebAssembly module. A change in a Vite config or build plugin, in a page, or in the threaded module can break the shipped web app or a threaded run without a failing check.

## Evidence

- The `CI` workflow runs one job per check group on every pull request and on pushes to `rust`, each through `dev/docker/run just check-group <group>` [.github/workflows/cli.yml#L65](../../.github/workflows/cli.yml#L65)
- The `clippy` group runs `lint-wasm` and `lint-wasm-mt`, and the `tests` group runs `test-wasm` [justfile#L124-L126](../../justfile#L124-L126). `test-wasm` builds the single-threaded module with `build-wasm "dev"` and runs the Bun smoke test `packages/app-wasm/smoke.ts` [justfile#L303-L304](../../justfile#L303-L304). `test-wasm-mt`, which runs the same smoke test under Node on the threaded module, is outside every gate [justfile#L309](../../justfile#L309)
- The `Web app` workflow runs `just build-web2` on every push to any branch and when a pull request opens or reopens, in the image of the WebAssembly target [.github/workflows/web.yml#L65](../../.github/workflows/web.yml#L65). It builds both modules in the dist profile, as shipped, and the Vite page, so a commit that breaks the production build of web2 gets a failing `Web app` check
- `build-web` builds the web app with `bun run build:web` [justfile#L677-L678](../../justfile#L677-L678). No check group and no workflow runs it. `test-ts` runs `bun test --parallel` and builds no app
- Vite loads each config in Node and leaves workspace imports of the config to Node. A config that Node cannot load breaks every build, dev server, and image of that app, while `just check` and `just test-ts` pass: the example-data plugin once imported the extensionless entry of `@neherlab/app-contracts`, and only a manual `just build-web2` showed it
- The threaded module depends on version pairs that must change together: the nightly of the threaded build and its `-Zbuild-std`, Emscripten and the threaded OpenBLAS archive (`binarylandia/build_openblas`), which links only with the Emscripten version that built it, and the pthread pool of the module. The production build of the `Web app` workflow links the threaded module but never runs it

## Fix direction

Add `build-web` and `test-wasm-mt` to a check group, so the production build of the web app and a run of the threaded module happen on every pull request.

## Validation

- A pull request that breaks the production build of the web app, or a run of the threaded module, fails a CI check

## Proposed implementation plan

1. Add `test-wasm-mt` to `checks_tests` in `justfile` (line 126), so the `tests` job of `.github/workflows/cli.yml` runs the threaded module under Node on every pull request. Confirm that the image of that job holds the nightly toolchain, Emscripten, and the threaded OpenBLAS archive that `build-wasm-mt` needs; if it does not, run `test-wasm-mt` in the `Web app` workflow (`.github/workflows/web.yml`), which already uses the image of the WebAssembly target.
2. Add `build-web` to `checks_typescript` (line 128), or as a step of the `Web app` workflow beside `just build-web2`, so the production build of `app-web` runs on every pull request.
3. `check_full` (line 131) concatenates the group variables, so `just check-all` then also runs both recipes without a further edit.
4. No unit test applies; this is CI configuration. Validate with `just check-group tests` and `just check-group typescript` on the branch, then see red: on a scratch branch, make a Vite config of `app-web` import an extensionless workspace entry and confirm that the CI job fails, then drop the branch.
