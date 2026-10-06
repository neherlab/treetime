# Grid point limit: time inference stops with an error above a configurable grid size

v1 limits the number of points of every probability grid that time inference creates. A grid operation that would need more points stops the run with an error that names the needed and allowed point counts. The default limit is 1 000 000 points. v0 has no limit.

## Problem

Time inference represents node-time and branch-time distributions on uniform grids. Some operations choose the grid spacing from their inputs, so the point count follows the ratio of a distribution's width to the finest input spacing, not the size of the tree:

- `convolution_function_function_fine()` in `packages/treetime-distribution/src/distribution_ops/convolve.rs` resamples both operands to the smaller of their two spacings before the FFT convolution
- `resample_to_mass_window()` in `packages/treetime-distribution/src/distribution_ops/mass_domain.rs` keeps the input spacing when it is finer than the window width divided by the 300 inference grid points
- multiplication and division in `multiply.rs` and `divide.rs` sample the support intersection at the finer operand spacing

A fixed clock rate far from the rate the data supports makes branch-time distributions very narrow. On `data/mpox/clade-ii/20` with `--clock-rate=0.5`, about 5000 times a realistic mpox rate, a branch-time distribution spans `6.5e-5` years with a spacing of `2.2e-7` years, and the node distribution it is convolved with spans 3.4 years. Resampling the node distribution to the branch spacing needs 15 583 861 points. Before the limit, that run did not finish within 240 s, the same run on `data/mpox/clade-ii/1000` was terminated by the operating system (SIGKILL), and multiplication and division silently switched to a coarser grid at 1 000 000 points. The web app shows a terminated run as a lost back end, with no message.

At realistic clock rates the largest grid of a run stays far below the limit:

| Dataset                                                    | Clock rate | Largest grid (points) |
| ---------------------------------------------------------- | ---------- | --------------------- |
| `flu/h3n2/20`, `ebola/20`, `zika/20`, `tb/20`, `ebola/362` | estimated  | below 10 000          |
| `mpox/clade-ii/20`, `mpox/clade-ii/100`                    | `1e-4`     | below 10 000          |
| `rsv/a/20`                                                 | estimated  | 12 385                |
| `flu/h3n2/200`                                             | estimated  | 11 075                |
| `mpox/clade-ii/1000`                                       | `1e-4`     | 115 156               |
| `mpox/clade-ii/1000`                                       | `5e-4`     | 103 072               |

The coarser-grid fallback of multiplication never fired on these runs.

## What v0 does

v0 (`packages/legacy/treetime`) has no grid point limit. On `data/mpox/clade-ii/20` with `--clock-rate 0.5 --max-iter 1`, v0's default joint time inference finishes in 19 s with 219 MB, and v0 with `--time-marginal always`, the marginal inference that v1 uses ([timetree-marginal-only-time-inference.md](timetree-marginal-only-time-inference.md)), fails after 23 s with `Unexpected behavior detected in multiply function when determining peak of function with y-values '[]'`. Stopping with an error matches the outcome of the reference algorithm on this input.

## What v1 does

- `struct MaxGridPoints` in `packages/treetime-grid/src/max_grid_points.rs` owns the limit and its validation: at least 1000 points (above the 300 points of the inference grids, which every run creates), default 1 000 000. The CLI flag, the settings file, the settings route, the server option, and the generated schemas reject the same values
- `MaxGridPoints::point_count()` compares the needed point count, computed as a floating-point number, with the limit before any integer conversion or allocation, and returns `GridPointLimitExceeded` above it. A count that overflows `usize` or is infinite, as when the spacing underflows, is reported as "more points than the limit"
- Every grid operation of time inference checks its count: the grid constructors `Grid::from_range_dx()`, `GridFn::resample_range_dx_clamped()`, and `GridFn::from_arrays_nonuniform()`; both resampled convolution operands and the convolution result of `len_a + len_b - 1` points; multiplication, products, and division; and the mass window of `rewindow_to_mass()` and `convolve_across_edge()`
- Multiplication and division stop with the same error at the limit instead of switching to a coarser grid, so no operation lowers the resolution without telling the user
- `run_timetree()` in `packages/treetime/src/timetree/inference/runner.rs` adds an explanation to the error: with a fixed `--clock-rate`, it names that rate as the likely cause; otherwise it names the estimated rate and the input dates. Both suggest raising `--max-grid-points` when enough memory is available
- The limit bounds single grids. It does not bound the memory of a run, which also grows with the number of nodes and the number of grids alive at the same time

## Where the limit comes from

- A run: `--max-grid-points` of `treetime timetree`, or `max_grid_points` in its config
- The desktop app: `analysis.max_grid_points` in the settings file of the app folder, which `PUT /api/app-settings/analysis` replaces. It fills runs whose config does not set the value
- A web server: `treetime-server --max-grid-points`. It fills runs whose config does not set the value and rejects a higher value, in the code preview and at run start. The default is the built-in default
- Otherwise: 1 000 000

`RunLimits::apply()` in `packages/app-commands/src/run_limits.rs` applies the app setting or the server option to every command whose config schema has a `max_grid_points` setting. `AppService` applies it to the config of the code preview before it builds the command line, and to the config of a run at its start, after input path confinement. The run record keeps the filled value, so a run reproduces from its command line and YAML, and "Edit and run again" shows the value as a changed setting.

## Public server value

`dev/deploy/hetzner/compose.yaml` runs `treetime-server` with two processing threads in a container limited to 3 GB. The server value is the largest of 1 000 000, 500 000, and 200 000 points for which every run of the following measurement stays below 1.5 GB, half of the container memory: `data/mpox/clade-ii/1000` with `-j 2 --max-iter=1`, at fixed clock rates from `1e-4` to `2e-3`, peak resident memory from `VmHWM` of the release build.

| Limit (points) | Clock rate | Outcome | Peak memory (MB) | Time (s) |
| --- | --- | --- | --- | --- |
| 1 000 000 | `1e-4` | finished | 864 | 55 |
| 1 000 000 | `2e-4` | finished | 546 | 18 |
| 1 000 000 | `5e-4` | finished | 686 | 40 |
| 1 000 000 | `7e-4` | limit error | 437 | 8 |
| 1 000 000 | `1e-3` | limit error | 435 | 7 |
| 1 000 000 | `2e-3` | limit error | 435 | 8 |
| 500 000 | `1e-4` | finished | 869 | 48 |
| 500 000 | `2e-4` | finished | 543 | 18 |
| 500 000 | `5e-4` | finished | 687 | 39 |
| 500 000 | `7e-4` | limit error | 436 | 9 |
| 500 000 | `1e-3` | limit error | 435 | 8 |
| 500 000 | `2e-3` | limit error | 434 | 9 |
| 200 000 | `1e-4` | finished | 867 | 48 |
| 200 000 | `2e-4` | finished | 540 | 19 |
| 200 000 | `5e-4` | finished | 686 | 40 |
| 200 000 | `7e-4` | limit error | 435 | 7 |
| 200 000 | `1e-3` | limit error | 435 | 8 |
| 200 000 | `2e-3` | limit error | 433 | 9 |

Every run stays below 0.9 GB at every candidate, so the server uses the largest candidate, 1 000 000. The limit does not change the peak memory in this measurement: runs at rates up to `5e-4` finish with grids far below 200 000 points, and runs from `7e-4` stop within 10 s on node distributions that reach millions to trillions of years into the past ([kb/issues/M-timetree-mass-window-grid-explodes-at-moderate-clock-rate.md](../issues/M-timetree-mass-window-grid-explodes-at-moderate-clock-rate.md)). Times are from a shared development machine and are approximate.

## Consequences

- A run that needs a grid above the limit ends with an actionable error instead of being terminated by the operating system
- A run that today used the coarser-grid fallback of multiplication or division now stops; no measured run reached the fallback
- The limit does not make every input succeed: the slowdown at spacing ratios below the limit is [kb/issues/M-timetree-marginal-dense-mpox-slow.md](../issues/M-timetree-marginal-dense-mpox-slow.md), and the mass windows that grow beyond any limit at moderate clock rates are [kb/issues/M-timetree-mass-window-grid-explodes-at-moderate-clock-rate.md](../issues/M-timetree-mass-window-grid-explodes-at-moderate-clock-rate.md)
