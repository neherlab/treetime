# Input-mode branch likelihoods: exact first cell, grid reach, and interpolation bounds

**Type**: Intentional divergence from v0.

With `--branch-length-mode=input`, the branch-length likelihood of an edge comes from its input branch length alone (v0 `BranchLenInterpolator`, the `branch_length_mode == 'input'` branch, [packages/legacy/treetime/treetime/branch_len_interpolator.py#L64-L102](../../packages/legacy/treetime/treetime/branch_len_interpolator.py#L64-L102)). v1 evaluates the same cost as v0 but represents it on its own grid (`fn InputBranchLikelihood::distribution()` [packages/treetime/src/timetree/inference/input_branch_length.rs](../../packages/treetime/src/timetree/inference/input_branch_length.rs)). This document records how that representation differs from v0 and which differences are accepted.

Approval: the user delegated the choice to the implementer, with the instruction to take the most scientifically sound and practical option in the direction of the design of the v0 author.

## Background

$k$ is the number of mutations on the branch (input length times the sequence length $L$), $\mu = 1/L$ the length of one mutation, $x$ the candidate branch length, and $c = 10^{-3}\mu$ the v0 offset `MIN_BRANCH_LENGTH * one_mutation`. v0 uses two costs:

- **Poisson** (input length below 0.05 substitutions per site): $y(x) = [x - k\mu \ln(x + c)]/\mu$
- **Gaussian** (input length of at least 0.05): $y(x) = \min\left(\tfrac{(\ell - x)^2}{2\sigma^2},\ 100 + \tfrac{|\ell - x|}{\sigma_t}\right)$, quadratic near the input length $\ell$ and linear beyond 100 nats

v0 marks the input mode as an approximation (the code comment "APPROXIMATE HANDLING OF BRANCH LENGTH PROPAGATOR"). Its offset $c$ only keeps the logarithm finite at $x = 0$, and the Gaussian cost turns linear at 100 nats to limit the penalty of a branch that conflicts with the clock. v0 evaluates the cost on a non-uniform grid of 125 points, dense near $x = 0$, up to `MAX_BRANCH_LENGTH = 4`, and interpolates between the points.

v1 represents time distributions on uniform grids ([distribution-tails-and-arithmetic.md](distribution-tails-and-arithmetic.md)). The input-mode grid starts at $x = 0$ with a step of a tenth of the peak width: $\sigma/10$ for the Gaussian and $\sqrt{\max(k\mu, \mu)\,\mu}/10$ for the Poisson cost. Between grid points the negative log likelihood is linear, except in cells that carry a law.

## Decisions

### First Poisson cell follows the exact cost

The term $k \ln(x + c)$ curves on the scale $c$, far below the grid step (about $100c$ for $k \le 1$). A straight line across the first cell $[0, h]$ differed from v0 by up to 30 nats ($k = 10$, $L = 1000$). A uniform grid at the step $c/4$ would need 80 000 to 3 800 000 points per branch.

The first cell of a Poisson branch therefore follows the exact cost through a cell law of `treetime-grid`:

- `struct EdgeCellLaw` [packages/treetime-grid/src/edge_cell_law.rs](../../packages/treetime-grid/src/edge_cell_law.rs) stores the grid edge $t_e$, the offset of the origin of a power law beyond the edge, and the exponent $b$. The input mode sets $t_e = 0$, the offset $c/r$ and $b = k$, with $r$ the branch rate that converts branch lengths into times. The offset, unlike an absolute origin, keeps the distance to the origin exact when a leaf date shifts the cell to absolute times
- `BoundaryBehavior::HardEdgeCell(EdgeCellLaw)` puts the law on the side of the grid edge. The support stays hard at the edge: no mass lies at negative branch lengths, as in v0 and in every other input-mode cell
- `struct EdgeCell` [packages/treetime-grid/src/edge_cell.rs](../../packages/treetime-grid/src/edge_cell.rs) evaluates the cell. With $d$ the distance from the origin and $d_a$ that of the inner cell end (the anchor), the cost is $y(d) = y_a - b \ln(d/d_a) + s\,(d - d_a)$. The slope $s$ follows from the two stored cell ends, so the law reproduces both ends exactly and carries the linear part $x/\mu$ of the cost. For the Poisson cost, $s$ equals $r/\mu$ and the law equals v0's cost everywhere in the cell. For $0 < k < 1$ in units of $\mu$ ($b < 1$ mutation), the cost has its minimum inside the cell, at $x = k\mu - c$
- The mass of a part of the cell is integrated by Gauss-Legendre panels from the densest point of the part outward, in the logarithm of the distance; the density is log-concave, so the integral stops 40 nats below the peak with an error below $e^{-40}$, and no sums cancel ([kb/algo/distribution.md](../algo/distribution.md#edge-cell-with-an-offset-power-law))
- The mass profile, the trimming to the mass window, the convolutions with points and ranges, negation, time scaling and the posterior statistics read the law. The maximum posterior region grows from a mode inside the cell. The shift by the exact date of a leaf keeps the law, so the messages of such leaves carry it into the products at their parents. An operation that moves the grid edge away from $t_e$, and a convolution of two functions ([kb/issues/N-distribution-fft-convolution-ignores-edge-cell-law.md](../issues/N-distribution-fft-convolution-ignores-edge-cell-law.md)), turn the side into a plain hard edge, and the cell becomes linear
- Products and quotients compose the laws of one side. At one edge the exponents add, and the offset is the exponent-weighted geometric mean of the offsets; the composed law is exact when the offsets agree and errs by at most $\sum_j b_j \ln^2(c_j/c^*)/8$ nats otherwise. Offsets differ when a relaxed clock gives two sibling leaves of the same date different branch rates: a rate ratio of 2 gives at most $0.03\,b$ nats for two factors of exponent $b$. A quotient inverts a product exactly. At different edges the law at the edge of the result stays, and the other factor is linear in the result's first cell with the interpolation error $b\,h^2/(8(\delta + c)^2)$, $\delta$ the distance between the edges: the bound of a later cell. An edge cell law and a hard-approach tail at the same point fail with an internal error: marginal-mode branch distributions carry hard-approach tails, input-mode ones edge cells, and one run uses one mode

The existing hard-cell law (`HardApproachLaw`) was not used with a grid start at $-c$: zero density at $-c$ puts a mass of up to $9.4 \cdot 10^{-4}$ of the branch below 0 ($k = 0$). Dropping the offset (the exact Poisson likelihood, zero density at 0) differs from v0 below about $20c$, and a branch that the dates force to length 0 would have no density.

### Grid reach of the Poisson cost

The Poisson grid ends where the cost is 100 nats above its minimum (`POISSON_GRID_NEG_LOG_SPAN`) and continues with a linear tail of the slope there, as v0's Gaussian cost turns linear at 100 nats. An exact tail to v0's `MAX_BRANCH_LENGTH = 4` would need a grid to 4 at the same step: 40 000 points ($k = 1$, $L = 1000$) to 220 000 points ($L = 30000$) per branch.

### Interpolation error of the Gaussian cost

Linear interpolation of a quadratic with curvature $1/\sigma^2$ errs by at most $h^2/(8\sigma^2)$ in a cell of width $h$. At the step $\sigma/10$ the bound is $1/800 = 1.25 \cdot 10^{-3}$ nats. This error is accepted. The grid ends where v0's cost turns linear, and the tail beyond it has v0's slope, so the representation is exact beyond the grid.

## Remaining differences from v0

Measured at the 134 grid points of the v0 interpolator (fixture `gm_input_branch_length_outputs.json` and the same grid for $k = 3$), within 100 nats of the peak, with both sides shifted to agree at $x = 0$. "Straight first cell" is the representation before the edge cell law. The largest difference with the edge cell lies in cell 1 for every case below $k = 499$; inside the first cell it is below $10^{-13}$ nats.

| Case | Straight first cell | Edge cell law | Bound in cell 1 ($k/8$) |
| --- | --- | --- | --- |
| Poisson $k = 0.05$, $L = 1000$ | 0.106 | $2.9 \cdot 10^{-3}$ | $6.3 \cdot 10^{-3}$ |
| Poisson $k = 1$, $L = 1000$ | 2.12 | 0.059 | 0.125 |
| Poisson $k = 3$, $L = 1000$ | 7.63 | 0.138 | 0.375 |
| Poisson $k = 10$, $L = 1000$ | 30.2 | 0.587 | 1.25 |
| Poisson $k = 30$, $L = 30000$ | 15.5 | 1.76 | 3.75 |
| Poisson $k = 499$, $L = 10000$ | $3.3 \cdot 10^{-3}$ | $3.3 \cdot 10^{-3}$ (cell 136) | |
| Gaussian (three fixture cases) | $9.7 \cdot 10^{-4}$ | $9.7 \cdot 10^{-4}$ | $1.25 \cdot 10^{-3}$ ($h^2/8\sigma^2$) |

Two differences remain, both recorded in [kb/issues/M-timetree-input-branch-likelihood-grid-diverges-from-v0.md](../issues/M-timetree-input-branch-likelihood-grid-diverges-from-v0.md):

- **Poisson cells after the first**: the curvature of $k \ln(x + c)$ in cell $i \ge 1$ is at most $k/(ih)^2$, so linear interpolation errs by at most $k/(8i^2)$ nats, independent of the step $h$. A finer uniform step does not reduce it; a grid that is non-uniform near 0, as in v0, does
- **Poisson cost beyond the grid**: the linear tail underestimates v0's convex cost by $k(r - 1 - \ln r)$, with $r$ the branch length divided by the grid end. At v0's last grid point ($x \approx 4$) the difference is 0.04% ($k = 0.05$) to 54% ($k = 499$) of v0's cost, which is more than 3900 nats above the peak in every case

## Effect on timetree results

The law changes the grid values of no input-mode distribution: the grid points hold v0's cost before and after it. It changes only what the operations read between the grid points of the first cell. In the time inference, the law reaches results through three channels:

- **Convolution window**: `fn convolution_output_window()` in [packages/treetime-distribution/src/distribution_ops/edge_convolution.rs](../../packages/treetime-distribution/src/distribution_ops/edge_convolution.rs) bounds each backward and forward message by the mass-bounded domain of its operands. The mass of the branch operand includes the exact first cell, so the quantile at the window edge moves, and with it the resampling grid of the message
- **Messages of leaves with exact dates**: the backward message of such a leaf is its branch distribution shifted to the leaf date, and it keeps the law. The products at the parent, the parent's posterior and its statistics read the law, also a mode inside the first cell
- **Range-function convolution**: the backward message of a leaf with an uncertain date (`2022-09-XX`, a date range) holds the mass of the branch distribution between two offsets at each output point (`fn convolution_range_function()` in [packages/treetime-distribution/src/distribution_ops/convolve.rs](../../packages/treetime-distribution/src/distribution_ops/convolve.rs)), and these masses include the exact first cell

Function-function convolution reads the edge cell by its grid values only and ends its result with a plain `Hard` side ([kb/issues/N-distribution-fft-convolution-ignores-edge-cell-law.md](../issues/N-distribution-fft-convolution-ignores-edge-cell-law.md)).

Measured with `treetime timetree --branch-length-mode=input --seed=1` (dev build, `--confidence` except for `data/rsv/a/20` and `data/sc2/4500`), before and after the law. "Shift" is the largest change of a node date; "of interval" divides each node's shift by the width of its 90% interval and reports the largest ratio:

| Dataset | Leaf dates | Changed nodes | Shift | Of interval | Log marginal likelihood of node times |
| --- | --- | --- | --- | --- | --- |
| `data/flu/h3n2/20` | exact | 18 of 37 | 4.68 days | 0.013 | -14.1617 to -14.1606 |
| `data/ebola/20` | exact | 19 of 39 | 2.22 days | 0.013 | -48.9739 to -49.0257 |
| `data/mpox/clade-ii/20` | 8 month or year ranges | 19 of 39 | 1.25 days | 0.023 | -69.0920 to -69.0281 |
| `data/mpox/clade-ii/100` | month or year ranges | 101 of 199 | 1.94 days | 0.049 | -415.1942 to -415.1923 |
| `data/rsv/a/20` | month ranges | 0 of 38 | 0 | 0 | -18.9681, unchanged |
| `data/sc2/4500` | exact | 4416 of 9561 (46 by more than 1 day, mean 0.08 days) | 1.66 days | not measured | -24510.9914 to -24511.4726 |

- **Clock rate**: unchanged in every case
- **Marginal mode**: outputs identical byte for byte before the law reached leaf messages (`data/flu/h3n2/20`, `data/ebola/20`); the golden masters of the marginal runner stay unchanged with it
- **Channels separated**: a build that ignores the law only in the convolution window gave outputs identical byte for byte to the run before the law on `data/flu/h3n2/20` and `data/ebola/20` while leaf messages ended with a plain hard edge. Its shifts, up to 1.3% of the interval, come from the placement of the resampling grid, not from the model ([kb/issues/M-timetree-node-dates-move-when-a-time-grid-gains-one-point.md](../issues/M-timetree-node-dates-move-when-a-time-grid-gains-one-point.md)). The leaf messages then moved dates on these datasets by less than 0.01 days, on `data/mpox/clade-ii/20` and `data/mpox/clade-ii/100` by up to 1.58 days (4.8% of the interval), and on `data/sc2/4500` by up to 1.35 days; the log marginal likelihood rose by 0.07 nats on `data/mpox/clade-ii/20` and by 1.8 nats on `data/sc2/4500`. Short branches with fewer than one mutation have the minimum of their cost inside the first cell, where a straight cell misses mass
- **`data/lassa/L/20`** fails in input mode with and without the law at the grid point limit ([kb/issues/M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md](../issues/M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md))

## Tests

- `packages/treetime-grid/src/__tests__/edge_cell.rs`: the cell law against the Poisson cost inside the first cell, its masses against quadrature of the cost and against closed forms of steep cells, its mode, the inverses of mass and value on both branches of the mode, the composition of laws, and the hard support beyond the edge
- `packages/treetime-distribution/src/distribution_ops/__tests__/test_edge_cell.rs`: the mass profile, negation, time scaling, trimming, convolutions, products and quotients with an edge cell, between grid points
- `packages/treetime/src/timetree/__tests__/test_posterior_edge_cell.rs`: the maximum posterior region and the survival function inside an edge cell, also around a mode inside the cell on both sides
- `packages/treetime/src/timetree/inference/__tests__/test_gm_input_branch_length.rs`: the first Poisson cell against v0 at the v0 grid points, the bound $k/(8i^2)$ in the later Poisson cells, the bound $h^2/(8\sigma^2)$ of the Gaussian cost, and its convergence at second order
