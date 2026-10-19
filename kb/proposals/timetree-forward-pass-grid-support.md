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

## Experiment results (2026-10-19)

Prototype on branch `claude-grid-trim` (base `26c1c8533`, uncommitted, switched by environment variables; patch in the
session scratchpad as `grid_trim_experiment.diff`). Trimming keeps a trimmed hard side `Hard` at the new edge (option
c); option (a), a `Linear` tail fitted at the new edge, put mass past child tip dates on the simulated data (a parent's
interval reached 0.15 years past its child's sampling date) and was dropped. Release build, `--jobs 1`, idle machine.

### Where to trim: the forward posterior only

| variant | what is trimmed by mass, hard sides included | result |
| --- | --- | --- |
| c | forward convolution output and forward posterior | `data/mpox/clade-ii/1000`: 152 internal nodes without a date (empty posterior), tips with uncertain dates moved by up to 6 years; `data/sc2/4500`: some intervals at 3% of base |
| p | forward posterior only | consistent with base on all datasets (below); without the spacing cap, `mpox` exceeds `--max-grid-points` |
| p + coarse | as p, convolution output spacing at most the coarser operand's | below |

Trimming the convolution output by its own mass is unsafe: the subtree distribution it is multiplied with can put
the product's mass exactly where the convolution was trimmed (the multiplicative hazard of the [empty
invariant](../decisions/distribution-tails-and-arithmetic.md#empty-invariant)). The backward messages must not be
trimmed for the same reason (a first prototype that trimmed every `fn convolve_across_edge()` call cut a subtree
message at 2024.22 where the posterior had mass up to 2024.38). Trimming the posterior is safe in exact arithmetic: a
parent's posterior mass removed below its $\epsilon$ quantile bounds the mass any child loses by $\epsilon$. Children
still get small grids, because their operand, the parent's cavity, comes from the trimmed posterior.

Trimming and the spacing cap have to land together, not in the order of the [issues](#issues): a trimmed posterior is
gridded at 300 points over its mass window, its fine spacing passes into the next convolution with a wide branch
distribution, and the forward grids grow (sc2: median forward convolution 9 989 to 13 349 points, maximum 240 050;
mpox: a 1.5-million-point grid).

### Speed (p + coarse)

| dataset | wall time base | p + coarse | forward pass, one round | peak memory |
| --- | --- | --- | --- | --- |
| validation `const/rep00/L10000_mu0.01`, fixed clock, 3 passes | 36.3 s | 5.9 s | 8.7 s to 0.3 s | 398 to 123 MB |
| `data/sc2/4500` | 47.3 s | 14.9 s | 9.8 s to 1.05 s | 1070 to 695 MB |
| `data/mpox/clade-ii/1000`, fixed clock | 5.6 s | 5.3 s | 0.33 s to 0.24 s | 414 to 410 MB |
| `data/ebola/362` | 1.0 s | 0.9 s | 0.08 s to 0.06 s | 60 MB |

Median forward grid (convolution / posterior points): simulated 85 424 / 84 245 to 459 / 461; sc2 9 989 / 9 620 to
449 / 300. On this base mpox no longer has the 30 000-point grids measured above, so it gains little. Validation
suite (2880 runs): run time ratio per variant 0.50-0.72, total 1300 s to 516 s on 12 cores.

### Accuracy

- Validation suite, 2880 runs against base: no failures, no node dated before its parent, median change of date RMSE
  per variant below $2\cdot10^{-4}$ years. 33 runs change their date RMSE by more than 0.05 years (16 worse, 17
  better), all at $L=10^3$, $\mu=3\cdot10^{-4}$, where the root date is off by 7-18 years in base too.
- Real data against base (no truth): every node dated; largest date change 0.005 years (ebola), 0.027 (mpox), 0.009
  (sc2); interval widths median 0.993-0.997 of base. On sc2, 97 of 5099 internal intervals are narrower by more than
  10%, all on a ladder of nodes pinned at the earliest tip date, where base itself dates a child interval 0.034 years
  before its parent's (one violation above 0.01 years in base, none trimmed).
- Calibration gate (240 sets, `--coalescent-opt`, estimated rate, `--time-marginal only-final`): 90% coverage
  unchanged within 0.013 wherever $\mu \ge 10^{-3}$ or $L=10^4$. At $L=10^3$, $\mu=3\cdot10^{-4}$ coverage falls
  from 0.830 to 0.743, in all 30 runs, with interval widths down 20-35% and the date RMSE unchanged. At this signal
  base posteriors are gridded coarsely over their whole hard support (one node: 68 points over five years), and base
  reports a 90% interval of 1.2 years for a posterior whose 5-95% quantiles span 0.15 years; trimmed posteriors are
  resolved at 300 points over their mass. The lower coverage is the model's under-coverage at low signal, which the
  coarse base grids had masked; it needs its own fix (rate and prior at low signal), not wider grids.

### Branch distributions (issue 4)

With the forward pass fixed, building the branch distributions is half the run time of the simulated set (0.9 s of
1.9 s per round). Of that, 93% is evaluating the per-site likelihood at the 300 grid branch lengths
(`fn evaluate_log_lh_on_grid()`), 5-6% building the edge contributions, about 1% re-windowing. It scales with the
mutations: 0.81 s per round at $\mu=10^{-2}$, 0.11 s at $\mu=10^{-3}$ ($L=10^4$, 400 tips). Directions: fewer
evaluation points with interpolation of the smooth log-likelihood, or grouping sites with the same pattern.

### Operand tails (issue 3)

Prototype: before each convolution, `Linear` tails of both operands are resampled onto the grid for 12 nats of decay
at the operand's spacing (`HardApproach` tails, the branch distributions' side near zero length, left as they are).
Calibration gate as above, without trimming: at $L=10^3$, $\mu=3\cdot10^{-4}$ coverage 0.830 to 0.817 and date RMSE
0.494 to 0.586; at a fixed true rate on one set the interval widths stay within 1% (median 0.994) and dates within
0.14 years, while with an estimated rate the probes of the rate step change (one set: rate $1.9\cdot10^{-4}$ to
$2.3\cdot10^{-4}$, true $3\cdot10^{-4}$). The prototype also fails on 31 of 240 runs, all with $\mu \ge 10^{-3}$,
$L=10^4$ or $\mu=10^{-2}$ (`NaN` when re-peaking a backward message), so its high-signal numbers are biased. No sign
that the dropped tails limit the low-signal coverage; the coarse-grid effect above is the larger one. Not pursued.

### Recommendation

Land posterior trimming (option c at the posterior only) together with the spacing cap
`dx = min(window / (GRID_POINTS - 1), max(dx_a, dx_b))`, both in the forward pass only, as one change; record in a
decision that the backward messages and the forward convolution output keep their hard support, and why. Then take
issue 4.

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
