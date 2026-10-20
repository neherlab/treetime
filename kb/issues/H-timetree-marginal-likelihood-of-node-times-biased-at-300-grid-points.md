# The marginal likelihood of node times is biased at 300 grid points

The clock rate step reads the rate from the log marginal likelihood of node times at three probe rates (`fn step_marginal_clock_rate()` [packages/treetime/src/timetree/optimization/marginal_clock_rate.rs](../../packages/treetime/src/timetree/optimization/marginal_clock_rate.rs)). The backward messages that this likelihood integrates have `GRID_POINTS = 300` points over a mass window that drops `EPS` $= 5 \cdot 10^{-4}$ of each message's mass ([packages/treetime/src/timetree/inference/runner.rs#L31-L38](../../packages/treetime/src/timetree/inference/runner.rs#L31-L38)). On simulated data the likelihood at 300 points is several nats above its converged value, and the bias depends on the rate enough to move the estimated rate by 1% to 14%.

## Evidence

Three sets of the simulated validation subset (`const` coalescent), log marginal likelihood against the number of message grid points and the mass window `EPS`, once at the reference rate without a prior, and once as the first evaluation of the command with `--coalescent-opt` and an estimated rate:

| Set, evaluation | 300 / $5 \cdot 10^{-4}$ | 600 | 1200 | 2400 | 4800 / $5 \cdot 10^{-4}$ | 4800 / $10^{-6}$ | Gap 300 to 4800 / $10^{-6}$ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| $L = 10^4$, $\mu = 3 \cdot 10^{-4}$ rep 0, keep root, no prior | -699.87 | -703.36 | -705.34 | -706.35 | -706.78 | -706.61 | 6.74 |
| $L = 10^4$, $\mu = 3 \cdot 10^{-4}$ rep 1, keep root, no prior | -932.16 | -934.82 | -936.37 | -937.21 | -937.59 | -937.46 | 5.30 |
| $L = 10^3$, $\mu = 10^{-3}$ rep 1, reroot, no prior | -473.26 | -479.50 | -483.18 | -484.95 | -485.85 | -485.59 | 12.33 |
| $L = 10^4$, $\mu = 3 \cdot 10^{-4}$ rep 0, `keep-opt-est` | -3022.06 | -3028.60 | -3032.48 | -3034.29 | -3035.22 | -3034.93 | 12.87 |
| $L = 10^4$, $\mu = 3 \cdot 10^{-4}$ rep 1, `keep-opt-est` | -3232.15 | -3236.72 | -3239.32 | -3240.82 | -3241.49 | -3241.24 | 9.09 |
| $L = 10^3$, $\mu = 10^{-3}$ rep 1, `reroot-opt-est` | -2683.88 | -2701.35 | -2709.57 | -2714.20 | -2716.36 | -2715.71 | 31.83 |

- **Convergence order**: the gap halves with each doubling of the points, which is first order in the spacing. The trapezoid rule of smooth messages is second order, so a first-order source dominates; [M-distribution-fft-convolution-hard-edge-first-order-error.md](M-distribution-fft-convolution-hard-edge-first-order-error.md) describes one such source
- **Mass window**: `EPS` changes the value by less than the grid does (at 4800 points, $5 \cdot 10^{-4}$ against $10^{-6}$: 0.1 to 0.7 nats)
- **Before the forward-posterior trim**: the code before the trim and the mass model has the same gaps (6.68, 5.24, 12.33, 12.67, 8.96, 30.41) and differs from the current code by 0.06 to 1.4 nats at 300 points
- **Rate**: at 4800 points the two codes agree within 0.13% in rate and 0.017 years in root error. The rate step at 300 points ends 1% to 14% away from the rate at 2400 and 4800 points, on both codes (`keep-opt-est` $L = 10^4$, $\mu = 3 \cdot 10^{-4}$ rep 0: $2.691 \cdot 10^{-4}$ against $2.721 \cdot 10^{-4}$; `reroot-opt-est` $L = 10^3$, $\mu = 10^{-3}$ rep 1: $9.72 \cdot 10^{-4}$ against $1.133 \cdot 10^{-3}$, true rate $10^{-3}$)

## Impact

- The estimated clock rate under `--clock-rate` absent and a coalescent prior, and through it every node date of such runs
- Changes of the discretization that leave the 4800-point result unchanged move the mean root error of two cells of the validation subset (2 replicates each) across the accuracy limit $\max(1.1 \cdot \text{base}, \text{base} + 0.05)$ years: `keep-opt-est` $L = 10^4$, $\mu = 3 \cdot 10^{-4}$ from 0.528 to 0.615 years (limit 0.581), and `reroot-opt-est` $L = 10^3$, $\mu = 10^{-3}$ from 1.828 to 2.043 years (limit 2.011)
- The reported log marginal likelihood (`log_lh_time`)

## Fix direction

Find the first-order source by refining one part of the backward pass at a time (FFT spacing, message grid, the hard edges of branch distributions near zero length, the coalescent factor), then fix it there. Raising `GRID_POINTS` alone is expensive: the untrimmed pipeline at 4800 points takes 1851 s of CPU and 7.9 GB on `data/sc2/4500`, against 164 s and 1.0 GB at 300 points.

> [!IMPORTANT]
> **Investigation required.** Which operation contributes the first-order term. The tests of the marginal likelihood and of the clock rate step pass at $10^{-6}$ against closed forms on small fixtures with fine branch grids, so the source appears with the coarse 300-point messages of real runs.

## Validation

On the sets above, the log marginal likelihood at the default grid lies within 0.5 nats of the 4800-point value, and the rate after the first step within 1% of it.
