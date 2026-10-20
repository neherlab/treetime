# The dataset size that web2 supports is not documented

The WebAssembly back end of web2 runs in a 32-bit address space without threads. The module grows its memory up to 4 GiB (`-sMAXIMUM_MEMORY=4GB` in `.cargo/config.toml`), and a run computes on one thread. Neither the documentation nor the interface tells a user which datasets fit, or that the CLI is faster on several threads.

## Evidence

Measured with the release profile of `wasm32-unknown-emscripten` in Node 24, against the native release CLI on the same machine with the configuration that web2 resolved for the same run:

| Run                            | web2   | CLI, 1 thread | CLI, 16 threads | web2 memory | CLI peak RSS |
| ------------------------------ | ------ | ------------- | --------------- | ----------- | ------------ |
| `sc2/4500` ancestral           | 5.6 s  | 3.2 s         | 1.8 s           | 1215 MiB    | 875 MiB      |
| `sc2/4500` timetree            | 46.1 s | 42.8 s        | 6.7 s           | 1227 MiB    | 1054 MiB     |
| `mpox/clade-ii/2000` ancestral | 11.3 s | 6.2 s         | 3.5 s           | 1695 MiB    | 1700 MiB     |
| `mpox/clade-ii/2000` timetree  | 23.4 s | 20.8 s        | 4.7 s           | 1476 MiB    | 1062 MiB     |

- web2 memory is the size of the linear memory of the module after the run; it never shrinks, so it is the peak. The Node process that ran the module peaked at 1.8 to 2.7 GiB, because the in-memory file system of Emscripten keeps the input and output files in JavaScript memory outside the linear memory
- On one thread, web2 takes 1.1 times the time of the CLI for `timetree` and 1.8 times for `ancestral`; against the CLI on 16 threads it takes 3 to 7 times as long
- The run outputs of `mpox/clade-ii/2000` ancestral are 1.1 GB, among them `ancestral.augur-node-data.json` with 785 MB, all held in the in-memory file system
- Browsers: Chrome and Firefox allow 4 GiB of `wasm32` memory; Safari and iOS Safari have no Memory64, and iOS Safari allows much less memory per page in practice ([emscripten#19144](https://github.com/emscripten-core/emscripten/issues/19144))

> [!IMPORTANT]
> **Investigation required.** The measurements above ran in Node, one run per instance. Measure the memory of the worker in Chrome, Firefox, and Safari, and the memory after several large runs in one session, because the linear memory never shrinks and the run store keeps every output file in memory.

> [!IMPORTANT]
> **Investigation required.** web2 applies no grid point limit of its own: `browser_config()` in `packages/app-wasm/src/backend.rs` sets no `max_grid_points` and no app settings, so a run without `max_grid_points` uses the built-in 1 000 000 points ([kb/decisions/distribution-grid-point-limit.md](../decisions/distribution-grid-point-limit.md)). An allocation that fails in the module aborts it, and the abort loses every run of the page ([M-web2-wasm-panic-loses-the-back-end-state.md](M-web2-wasm-panic-loses-the-back-end-state.md)). Measure the peak linear memory of `mpox/clade-ii/1000` `timetree` with `--max-iter=1` at fixed clock rates from `1e-4` to `2e-3` in the module, as the decision measured the public server, and check whether any run reaches 4 GiB below the limit error. Only such a run justifies a lower default for web2. That default would fill runs without a value and accept higher values, unlike `RunLimits::Server`, which rejects them.

## Fix direction

Document the dataset sizes that web2 supports, and show the CLI as the tool for larger datasets and for faster runs on several threads.

## Validation

- The documentation of web2 states the supported dataset size, measured in the browsers that web2 supports
