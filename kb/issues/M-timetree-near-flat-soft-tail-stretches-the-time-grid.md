# A soft tail with a slope near zero stretches a time grid far beyond the data

A `Linear` tail law with neg-log slope $k$ holds the mass $p_\text{edge} / |k|$ beyond its grid edge, where $p_\text{edge}$ is the density at the edge, and its mass quantiles lie at distances of order $\ln(1/\epsilon) / |k|$ from the edge for a mass fraction $\epsilon$. When an operation produces a decaying slope close to zero, a mass window or a mass quantile of that distribution lies years to trillions of years away from the data, and the next grid at the operand spacing exceeds the grid point limit or takes most of the run time. No code bounds the slope of a soft tail from below.

## Observed cases

- **Leaf range message** (fixed): on `data/mpox/clade-ii/1000` with `--clock-rate=1e-3 --max-iter=1`, the range convolution that shifted a function to the midpoint of a date range left a leaf message with the left slope $-4.6 \cdot 10^{-12}$; its $\epsilon$ quantile lay at $-1.66 \cdot 10^{12}$ years, and the run stopped on a grid over \[-1.66e12, 2022.62\] that needed $2.1 \cdot 10^{16}$ points. The exact range convolution ([kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md#exact-point-and-range-convolution)) removed this source; the run completes at the rates $5 \cdot 10^{-4}$, $7 \cdot 10^{-4}$, $10^{-3}$ and $2 \cdot 10^{-3}$, with 35 nodes keeping their given date, as before at $5 \cdot 10^{-4}$
- **Forward convolution output**: on `data/flu/h3n2/200` at `--clock-rate=3.058e-3 --keep-root` with a trim at $5 \cdot 10^{-4}$, the forward convolution of `NODE_0000157` (parent of 19 points over \[1972.891, 1972.994\]) had a soft right tail so shallow that its `EPS` mass window reached \[1972.89, 3523.67\], 1550 years, at the branch spacing: a grid of 99 095 points and an FFT of 589 824 against 12 288 without the trim
- **Polytomy resolution**: `data/zika/86` with `--coalescent-opt --resolve-polytomies --seed 7 --max-iter 6` stopped in round 2 on a grid over \[0, 2949.41\] (1 062 601 points). An empty sparse observation of a polytomy merger node caused it; the run completes since `fn reconcile_topology()` receives the merges of the resolution (`PolytomyResolution.merges`, `packages/treetime/src/partition/marginal/sparse/partition.rs`). Whether a near-flat tail carried the grid to year 0 is not shown
- **Soft window clamp refuted**: clamping the soft side of a mass window to a multiple of the input grid width failed 4 tests and moved flu dates by up to 13.5 years, because the far mass is real tail mass that the FFT drops later

> [!IMPORTANT]
> **Decision required.** How to guard against a soft tail whose slope is near zero. Options:
>
> - **Lower bound on the slope of a soft tail law**, for example from the width of the distribution's grid: changes the mass of the tail, which the log marginal likelihood and the posterior statistics read, and needs a rule for the bound
> - **Window cap relative to the input grid**: a mass window ends at most a fixed multiple of the operand's width from its grid; the refuted clamp above is one such rule, so a cap needs a measurement that it changes no result where the tail mass is real
> - **Error that names the node**: stop with an error that names the node and the tail slope when a mass window exceeds a bound, instead of the grid point error with its advice about the clock rate
>
> Evidence: the cases above; the narrow-branch case of [M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md](M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md) is a different mechanism (a fine spacing, not a wide window).

## Locations

- `fn mass_bounded_domain()` [packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L310-L330](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L310-L330): the soft side of a mass window from the tail quantile
- `fn SoftTailLaw::fit()` [packages/treetime-grid/src/soft_tail_law.rs#L13](../../packages/treetime-grid/src/soft_tail_law.rs#L13): clamps a rising slope to flat, with no lower bound on a decaying one
- `fn reconstruct_neg_log_tails()` [packages/treetime-distribution/src/distribution_ops/convolve.rs#L355](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L355): the secant tail of an FFT result

## Validation

- A test with a soft tail of slope $10^{-12}$ that enters a convolution: the run either keeps the grid below the limit with the decided guard, or stops with the decided error
- `data/mpox/clade-ii/1000` at the four rates above and `data/flu/h3n2/200` at the fixed rate keep their largest grid below 100 000 points
