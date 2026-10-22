# Input-mode branch likelihoods diverge from v0 between grid points and beyond the grid

## Problem

With `--branch-length-mode=input`, `fn InputBranchLikelihood::distribution()` ([packages/treetime/src/timetree/inference/input_branch_length.rs](../../packages/treetime/src/timetree/inference/input_branch_length.rs)) samples v0's branch-length likelihood (`branch_len_interpolator.py`, the `branch_length_mode == 'input'` branch) on a uniform grid. Between grid points, the distribution interpolates the negative log likelihood linearly. Beyond the last grid point, it continues linearly with the slope of the likelihood there. At the grid nodes, the distribution equals the v0 formula. The golden-master tests in `packages/treetime/src/timetree/inference/__tests__/test_gm_input_branch_length.rs` assert this, together with the grid end and the tail slope of the Gaussian case.

v0 evaluates the formula on its own non-uniform grid, which is dense near branch length 0 and reaches `MAX_BRANCH_LENGTH = 4`. The differences below are measured at the v0 grid points of the fixtures `gm_input_branch_length_outputs.json`, with each side shifted to 0 at the v0 grid point of the minimum. $k$ is the number of mutations on the branch, $L$ the sequence length, and $c = 10^{-3}/L$ the v0 offset `MIN_BRANCH_LENGTH * one_mutation`. No test asserts these differences.

## Measurements

The grid step is a tenth of the peak width (`GRID_POINTS_PER_WIDTH`).

- **Gaussian** (branches of at least 0.05 substitutions per site): the largest difference is $9.7 \cdot 10^{-4}$ nats, on both tails and in v0's linear tail out to $x = 5.8$. This is within the bound of linear interpolation, $h^2/8 = 1.25 \cdot 10^{-3}$ nats at a step of $\sigma/10$
- **Poisson near 0**: the term $k \ln(x + c)$ curves on the scale $c$, far below the grid step (about $70c$ for $k = 1$, $L = 1000$). In cell $i \geq 1$, the bound of linear interpolation is about $k/(8 i^2)$ nats and does not shrink with a finer step. Largest differences at the v0 points within 100 nats of the peak:

  | Case | Linear first cell | Exact first cell |
  | --- | --- | --- |
  | $k = 0.05$, $L = 1000$ | 0.106 | 0.003 |
  | $k = 1$, $L = 1000$ | 2.12 | 0.059 |
  | $k = 3$, $L = 1000$ | 7.63 | 0.138 |
  | $k = 10$, $L = 1000$ | 30.2 | 0.587 |
  | $k = 30$, $L = 30000$ | 15.5 | 1.76 |
  | $k = 499$, $L = 10000$ | 0.003 | 0.003 |

  The code has the linear first cell. The column "exact first cell" was measured with a script that follows $(x + c)^k$ exactly in cell 0
- **Poisson away from 0, inside the grid**: at most $3.3 \cdot 10^{-3}$ nats
- **Poisson beyond the grid**: the grid ends 100 nats above the peak (`POISSON_GRID_NEG_LOG_SPAN`). The Poisson cost is convex, so the linear tail underestimates v0 by $k(r - 1 - \ln r)$, where $r$ is the branch length divided by the grid end. At v0's last grid point ($x \approx 4$) the relative difference is 0.8% ($k = 1$), 6.4% ($k = 10$), 16% ($k = 30$, $L = 30000$) and 54% ($k = 499$). All these points lie more than 100 nats above the peak

## Open questions

> [!IMPORTANT]
> **Decision required.** Representation near branch length 0. The options:
>
> - An exact first cell with a new cell law in `treetime-grid`: a power law whose origin lies at $-c$, below the grid start. The support stays hard at 0
> - An exact first cell with the existing hard-cell law: the grid starts at $-c$ with zero density there. The mass below 0 is $9.4 \cdot 10^{-4}$ of the branch for $k = 0$, $6.5 \cdot 10^{-4}$ for $k = 0.05$ and $4.7 \cdot 10^{-7}$ for $k = 1$
> - A non-uniform grid near 0, as in v0. It removes the errors in the cells after the first one, but it changes the uniform `GridFn` that convolution, products and the mass profile use
> - The exact Poisson likelihood without v0's offset $c$: zero density at 0, as in marginal mode. It differs from v0 below about $20c$, and a branch that the dates force to length 0 then has no density
>
> A dense uniform grid at step $c/4$ needs 80 000 to 3 800 000 points per branch.

> [!IMPORTANT]
> **Decision required.** Extent of the Poisson grid. An exact tail out to v0's `MAX_BRANCH_LENGTH = 4` needs a grid out to 4 at the same step: 40 000 points ($k = 1$, $L = 1000$) to 220 000 points ($L = 30000$) per branch. The grid now ends 100 nats above the peak.

## Validation

- At the v0 grid points of the fixtures, within 100 nats of the peak, every Poisson case differs from v0 by less than the Gaussian bound
- Beyond 100 nats, the difference follows whichever decision is made on the extent of the Poisson grid
