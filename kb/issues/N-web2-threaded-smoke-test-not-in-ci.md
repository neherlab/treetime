# The smoke test of the threaded WebAssembly module does not run in CI

The `Web app` workflow builds web2 for production on every push, which links the threaded WebAssembly module, but no check runs a computation on the threaded module. A change in the threaded module or in one of the versions it depends on can break every threaded run of web2 without a failing check.

## Evidence

- The `CI` workflow runs one job per check group on every pull request and on pushes to `rust`, each through `dev/docker/run just check-group <group>` [.github/workflows/cli.yml#L65](../../.github/workflows/cli.yml#L65)
- The `clippy` group runs `lint-wasm` and `lint-wasm-mt`, and the `tests` group runs `test-wasm`, which builds the single-threaded module with `build-wasm "dev"` and runs the Bun smoke test `packages/app-wasm/smoke.ts`. `test-wasm-mt`, which runs the same smoke test under Node on the threaded module with 2 and 8 threads, is outside every gate (`justfile`, recipes `test-wasm` and `test-wasm-mt`)
- The `Web app` workflow runs `just build-web2` on every push to any branch and when a pull request opens or reopens, in the image of the WebAssembly target [.github/workflows/web.yml#L65](../../.github/workflows/web.yml#L65). It builds both modules in the dist profile, as shipped, and the Vite page, so a commit that breaks the production build of web2, its Vite config, or a build plugin gets a failing `Web app` check
- The threaded module depends on version pairs that must change together: the nightly of the threaded build and its `-Zbuild-std`, Emscripten and the threaded OpenBLAS archive (`binarylandia/build_openblas`), which links only with the Emscripten version that built it, and the pthread pool of the module. The production build links the threaded module but never runs it

## Fix direction

Add `test-wasm-mt` to a check group, so a run of the threaded module happens on every pull request.

## Validation

- A pull request that breaks a run of the threaded module fails a CI check
