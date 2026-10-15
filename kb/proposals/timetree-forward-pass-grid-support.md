# Proposal: bound timetree grids by probability mass, not by inherited support

## Problem

Most of the timetree run time goes to grid points that carry no probability. In the forward pass every node posterior keeps a hard lower edge inherited from the root and is gridded at the finest spacing seen anywhere in the tree, so its grid is one to two orders of magnitude longer than the region holding its mass. Every forward convolution, division and multiplication then works on that grid.

## Measurements

Base `94f921f8e7`, release build, `--jobs 1`, one JSON record per distribution operation at the call sites in `backward_pass.rs`, `forward_pass.rs` and `branch_length_likelihood.rs`. Mass core: the interval holding all but $10^{-3}$ of the grid mass ($5 \cdot 10^{-4}$ per side). Datasets: `data/mpox/clade-ii/1000` (fixed clock `5.7e-5`), `data/sc2/4500`, `data/ebola/362`, and the validation suite's `const/rep00/L10000_mu0.01` (400 tips, fixed clock, `--keep-root`).

Share of distribution-operation time:

| Operation | mpox | simulated |
| --- | --- | --- |
| forward convolution | 62% | 61% |
| forward division (cavity) | 21% | 21% |
| forward multiplication | 12% | 15% |
| backward convolution | 5% | 3% |
| backward product, branch distributions | < 1% | < 1% |

Output grids, median over calls:

| Dataset | Forward grid (points) | Total / core points, forward ops | Total / core points, backward ops |
| --- | --- | --- | --- |
| mpox | 30 000 – 35 000 | 54 – 117 | 1.4 – 2.3 |
| simulated | 84 000 – 86 000 | 58 – 124 | 1.6 – 2.5 |
| sc2 | 9 000 – 9 500 | 17 – 40 | 1.4 – 2.2 |
| ebola | 1 000 – 1 500 | 3 – 6 | 1.3 – 2.3 |

The excess lies on the hard left side of the forward distributions, where the density at the grid edge is below $10^{-50}$ of the peak. Posterior lower edges sit a median 250–1200 standard deviations below the mean (34 on ebola).

Wall time of the simulated run: forward pass about 55%, building branch distributions from the alignment about 35%, backward pass and the rest about 10%.

## Mechanism

Checked on every forward call of all four datasets (exceptions: 3–6 calls per dataset):

1. Forward convolution: lower edge = parent cavity lower edge + branch lower edge (about 0). `fn convolve_across_edge()` clips the window to the result's own `x_min` on the hard side, and `fn lower_edge()` / `fn upper_edge()` in `mass_domain.rs` return the grid edge for `Hard` and `HardApproach` sides, so hard sides are never trimmed by mass.
2. Division keeps the dividend's lower edge.
3. Multiplication of the forward message (hard left) with the subtree distribution (soft left) spans the full support intersection at the finest operand spacing ([grid intersection contract](../decisions/distribution-tails-and-arithmetic.md#grid-intersection-contract)), so its lower edge is the forward message's.

Each posterior therefore reaches back to about the root's lower edge. Separately, spacing never coarsens: `fn resample_to_mass_window()` takes the smaller of the window spacing and the input spacing, and the input of a convolution is on the finer operand spacing.

## Projected effect

Forward FFT length (sum of operand lengths on the finer spacing), median / 90th percentile:

| Dataset | Now | Operands trimmed to their mass cores | Trimmed, on the coarser spacing |
| --- | --- | --- | --- |
| mpox | 34 676 / 48 482 | 719 / 1 126 | 452 / 652 |
| simulated | 85 626 / 92 223 | 1 971 / 3 335 | 417 / 956 |
| sc2 | 9 467 / 9 830 | 637 / 1 012 | 441 / 627 |
| ebola | 1 533 / 2 045 | 560 / 803 | 421 / 529 |

Maximum after trimming on the finer spacing: 52 732 (mpox), 26 636 (simulated).

## Design axes

### Where to trim

- **Convolution output window and forward posterior.** Trim both sides by mass in `fn convolve_across_edge()` and after the forward multiplication. Covers the measured cost; division then inherits trimmed inputs.
- **Every operation result.** Also trim multiplication and division outputs generally. Uniform, but changes backward results that are not a measured cost.

### What a trimmed hard side becomes

- **`Linear` tail fitted at the new edge.** Support becomes unbounded on that side; under the [empty invariant](../decisions/distribution-tails-and-arithmetic.md#empty-invariant) a soft side bridges gaps, so trimming cannot make two messages disjoint. Drops the hard bound, which at the measured edges has density below $10^{-50}$.
- **`Linear` tail truncated at the original bound.** Keeps the exact support. Needs a new `BoundaryBehavior` variant and its handling in multiplication, division, mass and tail composition.
- **`Hard` at the trimmed edge.** Moves a hard bound inward. A prototype that shifted a hard edge inward produced a node with an empty posterior on `data/sc2/4500` (disjoint supports).

### Output spacing

- **Finest operand (current).** Keeps resolution from narrow factors everywhere downstream.
- **Window / (`GRID_POINTS` − 1), capped at the coarser operand spacing.** Resolution follows the result's own width.

### Operand tails in the convolution

The FFT convolves grid values only; operand tails beyond the grid are dropped and the result's soft tail is refitted. For a branch-length distribution windowed at $\epsilon = 5 \cdot 10^{-4}$ the dropped exponential tail carries about 3% of the variance of a zero-mutation branch (0.4% of its mean).

- **Extend soft tails onto the grid** for a fixed number of decay lengths before the FFT.
- **Exponential-kernel recursion:** convolution with $p_e e^{-k t}$ against a piecewise-linear function is a first-order recursion, exact for the interpolant and $O(n)$ at any $k \cdot dx$.
- **Linearity split:** convolve a gridded core and analytic tail pieces separately and sum.

### Moment matching

Replacing an operand whose standard deviation is at most the other's spacing by a mean shift and the three-point kernel $[s/2, 1-s, s/2]$, $s = \sigma^2/dx^2$, is exact for mass, mean and variance. It applies to 0–6% of the measured convolutions, so it is not a performance lever while supports are untrimmed. Prototype with tail moments (`TailMoments`, `fn SoftTailLaw::moments()`, `fn HardApproachLaw::moments()`, `fn MassProfile::mean_variance()`, `fn MassProfile::mass_quantile()`) on the unmerged branch `moment-matched-convolution`.

## Issues

1. [M-timetree-forward-grids-keep-untrimmed-hard-support.md](../issues/M-timetree-forward-grids-keep-untrimmed-hard-support.md): trim hard sides by mass
2. [M-timetree-grid-spacing-inherits-finest-operand.md](../issues/M-timetree-grid-spacing-inherits-finest-operand.md): output spacing
3. [M-distribution-convolution-drops-operand-tails.md](../issues/M-distribution-convolution-drops-operand-tails.md): operand tails in the convolution
4. [M-timetree-branch-distribution-build-cost.md](../issues/M-timetree-branch-distribution-build-cost.md): the next bottleneck after 1 and 2

Order: 1, then 2, each measured on its own; 3 after its investigation; 4 after 1 and 2.

## Related

- [M-timetree-mass-window-grid-explodes-at-moderate-clock-rate.md](../issues/M-timetree-mass-window-grid-explodes-at-moderate-clock-rate.md): support growth in the backward pass at high clock rates
- [H-timetree-mass-sizing-node-times-break-downstream-invariants.md](../issues/H-timetree-mass-sizing-node-times-break-downstream-invariants.md): re-gridding shifted committed node times and broke pipeline tests
- [exponential-branch-length-tails.md](exponential-branch-length-tails.md)
