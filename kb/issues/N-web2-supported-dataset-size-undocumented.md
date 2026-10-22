# The dataset size that web2 supports is not documented

Each run of web2 computes in a job worker of its own, with its own instance of the WebAssembly module in a 32-bit address space ([kb/decisions/app-web2-wasm-worker-back-end.md](../decisions/app-web2-wasm-worker-back-end.md)). Each instance grows its memory up to 4 GiB (`-sMAXIMUM_MEMORY=4GB` in [.cargo/config.toml#L73](../../.cargo/config.toml#L73)), so each run has its own 4 GiB. Runs that compute at the same time share the memory that the browser gives the tab, and when they exceed it, the browser closes the tab with every run. The back-end worker holds the uploads and the outputs of every run of the page in its own instance. Neither the documentation nor the interface tells a user which datasets fit, how many large runs can compute at the same time, or when the CLI is faster.

## Evidence

Measured with the release profile of the single-threaded module in Node 24, before runs moved into job workers, against the native release CLI on the same machine with the configuration that web2 resolved for the same run:

| Run                            | web2   | CLI, 1 thread | CLI, 16 threads | web2 memory | CLI peak RSS |
| ------------------------------ | ------ | ------------- | --------------- | ----------- | ------------ |
| `sc2/4500` ancestral           | 5.6 s  | 3.2 s         | 1.8 s           | 1215 MiB    | 875 MiB      |
| `sc2/4500` timetree            | 46.1 s | 42.8 s        | 6.7 s           | 1227 MiB    | 1054 MiB     |
| `mpox/clade-ii/2000` ancestral | 11.3 s | 6.2 s         | 3.5 s           | 1695 MiB    | 1700 MiB     |
| `mpox/clade-ii/2000` timetree  | 23.4 s | 20.8 s        | 4.7 s           | 1476 MiB    | 1062 MiB     |

- web2 memory is the size of the linear memory of the module after the run; it never shrinks, so it is the peak. The Node process that ran the module peaked at 1.8 to 2.7 GiB, because the in-memory file system of Emscripten keeps the input and output files in JavaScript memory outside the linear memory
- On one thread, web2 took 1.1 times the time of the CLI for `timetree` and 1.8 times for `ancestral`; against the CLI on 16 threads it took 3 to 7 times as long
- The run outputs of `mpox/clade-ii/2000` ancestral are 1.1 GB, among them `ancestral.node-data.json` with 785 MB. A job worker holds them in its file system until it has sent them, and the back-end worker holds them for the rest of the page session
- A trial of the threaded module in Chrome measured `sc2/4500` timetree at 15.8 s with threads against 1 min 44 s single-threaded

Times and peak linear memory of the job module by the URL parameter `jobs` (`jobs=1` is the single-threaded module), measured with the dist builds under Node, compute time only, in a container limited to 20 CPUs on a machine that other jobs loaded with 20 to 36 busy threads, so the times are noisy. On the same machine without load, the single-threaded module computed `sc2/4500` timetree in 18.7 s:

| `jobs` | `sc2/4500` timetree | `mpox/clade-ii/2000` ancestral |
| ------ | ------------------- | ------------------------------ |
| 1      | 45.5 s, 611 MiB     | 13.3 s, 1685 MiB               |
| 2      | 39.5 s, 611 MiB     | 12.7 s, 1634 MiB               |
| 4      | 34.3 s, 611 MiB     | 16.4 s, 1708 MiB               |
| 8      | 24.1 s, 611 MiB     | 14.6 s, 1668 MiB               |
| 16     | 18.2 s, 611 MiB     | 14.7 s, 1681 MiB               |

- At 16 threads, `sc2/4500` timetree is 2.5 times as fast as on one thread, and `mpox/clade-ii/2000` ancestral is not faster. Threads do not change the peak linear memory of a run
- In the threaded module every allocation takes the global lock of dlmalloc, the default allocator of Emscripten. The hot paths still allocate on every call: one `expQt` matrix per call and `.t().to_owned()` copies in [packages/treetime/src/partition/marginal/sparse/message.rs#L145-L169](../../packages/treetime/src/partition/marginal/sparse/message.rs#L145-L169), two `dot` results and a `GTR` clone per edge in [packages/treetime/src/partition/optimize/dense.rs#L10-L11](../../packages/treetime/src/partition/optimize/dense.rs#L10-L11), `Vec`s per edge in `packages/treetime/src/partition/optimize/sparse.rs`, and the `map_collect` of `divide_out` in `packages/treetime/src/partition/marginal/shared/normalize.rs`

> [!IMPORTANT]
> **Investigation required.** Whether the allocations above, or another cause, limit the speed-up of threaded runs. Profile `mpox/clade-ii/2000` ancestral in the threaded module with 16 threads on a machine without other load.

- Browsers: Chrome and Firefox allow 4 GiB of `wasm32` memory; Safari and iOS Safari have no Memory64, and iOS Safari allows much less memory per page in practice ([emscripten#19144](https://github.com/emscripten-core/emscripten/issues/19144)). A browser that cannot create a shared memory of 4 GiB runs every job on the single-threaded module

> [!IMPORTANT]
> **Investigation required.** The measurements above ran in Node, one run per instance. Measure the memory of the tab in Chrome, Firefox, and Safari with one large run, with several large runs computing at the same time, and after several large runs in one session, because the back-end worker keeps every upload and output file of the page in memory.

> [!IMPORTANT]
> **Investigation required.** web2 applies no grid point limit of its own: `browser_config()` [packages/app-wasm/src/backend.rs#L131-L139](../../packages/app-wasm/src/backend.rs#L131-L139) sets no `max_grid_points` and no app settings, so a run without `max_grid_points` uses the built-in 1 000 000 points ([kb/decisions/distribution-grid-point-limit.md](../decisions/distribution-grid-point-limit.md)). An allocation that fails in a job worker aborts its module and fails that run with the allocation size as the error. Whether a run within the grid point limit can reach the 4 GiB memory limit is not known. Measure the peak linear memory of `mpox/clade-ii/1000` `timetree` with `--max-iter=1` at fixed clock rates from `1e-4` to `2e-3` in the module, as the decision measured the public server, and check whether any run reaches 4 GiB below the limit error. Only such a run justifies a lower default for web2. That default would fill runs without a value and accept higher values, unlike `RunLimits::Server`, which rejects them.

## Fix direction

Document the dataset sizes that web2 supports and the number of large runs that can compute at the same time, and show the CLI as the tool for larger datasets.

## Validation

- The documentation of web2 states the supported dataset size, measured in the browsers that web2 supports
