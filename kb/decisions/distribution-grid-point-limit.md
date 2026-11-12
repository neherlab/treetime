# Grid point limit: time inference stops with an error above a configurable grid size

v1 limits the number of points of every probability grid that time inference creates. A grid operation that would need more points stops the run with an error that names the needed and allowed point counts. The default limit is 1 000 000 points. v0 has no limit.

## Problem

Time inference represents node-time and branch-time distributions on uniform grids. Some operations choose the grid spacing from their inputs, so the point count follows the ratio of a distribution's width to the finest input spacing, not the size of the tree:

- `convolution_function_function_fine()` in `packages/treetime-distribution/src/distribution_ops/convolve.rs` resamples both operands to the smaller of their two spacings before the FFT convolution
- `resample_to_mass_window()` in `packages/treetime-distribution/src/distribution_ops/mass_domain.rs` keeps the input spacing of a backward message when it is finer than the window width divided by the 300 inference grid points; a forward convolution output takes at most the coarser operand spacing unless its posterior needs a finer one ([timetree-forward-posterior-trimming.md](timetree-forward-posterior-trimming.md#spacing-cap-of-forward-convolutions))
- multiplication and division in `multiply.rs` and `divide.rs` sample the support intersection at the finer operand spacing

A fixed clock rate far from the rate the data supports makes branch-time distributions very narrow. On `data/mpox/clade-ii/20` with `--clock-rate=0.5`, about 5000 times a realistic mpox rate, a branch-time distribution spans `6.5e-5` years with a spacing of `2.2e-7` years. The forward message keeps that spacing, and the product of the forward message with a subtree distribution over \[2019.14, 2022\] needed 13 175 607 points on the code of 2026-10-19 ([kb/issues/M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md](../issues/M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md)). Before the limit, that run did not finish within 240 s, the same run on `data/mpox/clade-ii/1000` was terminated by the operating system (SIGKILL), and multiplication and division silently switched to a coarser grid at 1 000 000 points. The web app shows a terminated run as a lost back end, with no message.

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

- `struct MaxGridPoints` in `packages/treetime-grid/src/max_grid_points.rs` owns the limit and its validation: at least 1000 points (above the 300 points of the inference grids, which every run creates) and at most 9 007 199 254 740 991 points ($2^{53} - 1$, the largest count whose indices all convert to `f64` exactly; `usize::MAX` on 32-bit targets), default 1 000 000. The CLI flag, the settings file, the settings route, and the generated schemas reject the same values
- `MaxGridPoints::point_count()` compares the needed point count, computed as a floating-point number, with the limit before any integer conversion or allocation, and returns `GridPointLimitExceeded` above it. A count that overflows `usize` or is infinite, as when the spacing underflows, is reported as "more points than the limit"
- Every grid operation of time inference checks its count: the grid constructors `Grid::from_range_dx()` and `GridFn::resample_range_dx_clamped()`; the exact range convolution; both resampled convolution operands and the convolution result of `len_a + len_b - 1` points; multiplication, products, and division; and the mass window of `rewindow_to_mass()` and `convolve_across_edge()`
- Multiplication and division stop with the same error at the limit instead of switching to a coarser grid, so no operation lowers the resolution without telling the user
- `run_timetree()` in `packages/treetime/src/timetree/inference/runner.rs` adds an explanation to the error. Each grid request of time inference records its operation, the node, and its operands (`GridRequest` in `packages/treetime/src/timetree/inference/grid_request.rs`), and the explanation names the operation and node, the span and spacing of the grid, and the operands that set them. It names the fixed or estimated clock rate as the likely cause only when a branch-time distribution sets the spacing, because the clock rate sets the width of those distributions and of no other operand. Every explanation suggests raising `--max-grid-points` when enough memory is available. The clock-rate probes of `step_marginal_clock_rate()` in `packages/treetime/src/timetree/optimization/marginal_clock_rate.rs`, which run the backward pass outside `run_timetree()` at the estimated rate and at rates above and below it, explain their errors with the same function: the explanation says that estimating the clock rate needed the grid and names the probe rate as the rate that may make branch-time distributions too narrow or too wide
- The limit bounds single grids. It does not bound the memory of a run, which also grows with the number of nodes and the number of grids alive at the same time

## Where the limit comes from

- A run: `--max-grid-points` of `treetime timetree`, or `max_grid_points` in its config
- The desktop app: `analysis.max_grid_points` in the settings file of the app folder, which `PUT /api/app-settings/analysis` replaces. It fills runs whose config does not set the value
- Otherwise: 1 000 000

`RunLimits::apply()` in `packages/app-backend/src/run_limits.rs` applies the app setting to every command whose config schema has a `max_grid_points` setting. `AppService` applies it to the config of the code preview before it builds the command line, and to the config of a run at its start, after input path confinement. The run record keeps the filled value, so a run reproduces from its command line and YAML, and "Edit and run again" shows the value as a changed setting.

## Consequences

- A run that needs a grid above the limit ends with an actionable error instead of being terminated by the operating system
- A run that today used the coarser-grid fallback of multiplication or division now stops; no measured run reached the fallback
- The limit does not make every input succeed: a narrow branch distribution at a clock rate far above the one the data supports still needs millions of points ([kb/issues/M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md](../issues/M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md))
