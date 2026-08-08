# Node-time convergence tolerance is not derived from the branch grid

`NODE_TIME_TOLERANCE_YEARS = 1e-2` in
[`packages/treetime/src/timetree/convergence/metrics.rs`](../../packages/treetime/src/timetree/convergence/metrics.rs)
is a fixed constant. Node times are quantized by the `BRANCH_GRID_SIZE = 300` branch-length grid,
whose step depends on the dataset, so a fixed tolerance is above the resolvable step on some
datasets and below it on others.

## Evidence

Per-round maximum node-time changes come out as integer multiples of a per-dataset quantum:

| dataset | quantum | tolerance relative to it |
| --- | --- | --- |
| `data/ebola/20` | about 0.0026 y | 1e-2 is about 4 steps |
| `data/mpox/clade-ii/1000` | about 0.016 y (about 5.8 days) | 1e-2 is below one step |

On mpox the loop can reach a bit-identical `log_lh_seq` while `max_time_change` keeps jittering at
0.016 / 0.032 / 0.048 / 0.080 — node times snapping between adjacent grid points. Runs on such
datasets exhaust `--max-iter` rather than reporting convergence, even though the answer has
stopped changing.

## Impact

Wasted iterations on datasets whose grid step exceeds the tolerance. No incorrect result: the
extra rounds do not move the answer, and `--max-iter` bounds the cost.

## Options

- Derive the tolerance from the grid step, e.g. one step or a small multiple, computed per run.
- Keep a fixed floor at the resolution the data carry (see
  [timetree-convergence-on-node-times.md](../decisions/timetree-convergence-on-node-times.md)) and
  take the maximum of the two.
- Switch the criterion to the objective, with time movement as a guard rather than the primary
  test. This depends on
  [M-timetree-convergence-metric-deficiencies.md](M-timetree-convergence-metric-deficiencies.md),
  since the reported total is not currently comparable across rounds.

## Related issues

- [M-timetree-branch-grid-uniform-resolution.md](M-timetree-branch-grid-uniform-resolution.md)
- [M-timetree-convergence-metric-deficiencies.md](M-timetree-convergence-metric-deficiencies.md)
