# Proposal: bound timetree grids by probability mass, not by inherited support

## Status

The recommendation below is implemented and recorded in [kb/decisions/timetree-forward-posterior-trimming.md](../decisions/timetree-forward-posterior-trimming.md): every refined forward posterior is cut at its $10^{-9}$ mass quantiles (`TRIM_EPS`), a cut hard side becomes `Hard` at the cut, and forward convolution outputs take at most the coarser operand spacing unless the posterior they form needs a finer one. Backward messages and forward convolution outputs keep their support. The cut slices the posterior's own grid and does not resample it. The measurements and the mechanism below describe the code before that change (base `94f921f8e7`).

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
| mpox | 30 000 - 35 000 | 54 - 117 | 1.4 - 2.3 |
| simulated | 84 000 - 86 000 | 58 - 124 | 1.6 - 2.5 |
| sc2 | 9 000 - 9 500 | 17 - 40 | 1.4 - 2.2 |
| ebola | 1 000 - 1 500 | 3 - 6 | 1.3 - 2.3 |

The excess lies on the hard left side of the forward distributions, where the density at the grid edge is below $10^{-50}$ of the peak. Posterior lower edges sit a median 250-1200 standard deviations below the mean (34 on ebola).

Wall time of the simulated run: forward pass about 55%, building branch distributions from the alignment about 35%, backward pass and the rest about 10%.

## Mechanism

Checked on every forward call of all four datasets (exceptions: 3-6 calls per dataset):

1. Forward convolution: lower edge = parent cavity lower edge + branch lower edge (about 0). `fn convolve_across_edge()` clips the window to the result's own `x_min` on the hard side, and the mass window returns the grid edge for `Hard` and `HardApproach` sides, so hard sides are never trimmed by mass.
2. Division keeps the dividend's lower edge.
3. Multiplication of the forward message (hard left) with the subtree distribution (soft left) spans the full support intersection at the finest operand spacing (grid intersection contract, [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md#grid-intersection-contract)), so its lower edge is the forward message's.

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
- **Every operation result.** Also trim multiplication and division outputs. Uniform, but changes backward results that are not a measured cost.

### What a trimmed hard side becomes

- **`Linear` tail fitted at the new edge.** Support becomes unbounded on that side; under the empty invariant ([kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md#empty-invariant)) a soft side bridges gaps, so trimming cannot make two messages disjoint. Drops the hard bound, which at the measured edges has density below $10^{-50}$.
- **`Linear` tail truncated at the original bound.** Keeps the exact support. Needs a new `BoundaryBehavior` variant and its handling in multiplication, division, mass and tail composition.
- **`Hard` at the trimmed edge.** Moves a hard bound inward. A prototype that shifted a hard edge inward produced a node with an empty posterior on `data/sc2/4500` (disjoint supports).

### Output spacing

- **Finest operand (current).** Keeps resolution from narrow factors everywhere downstream.
- **Window / (`GRID_POINTS` - 1), capped at the coarser operand spacing.** Resolution follows the result's own width.

### Operand tails in the convolution

The FFT in `fn convolution_function_function_fine()` convolves grid values only; operand tails beyond the grid are dropped and the result's soft tail is refitted. For a branch-length distribution windowed at $\epsilon = 5 \cdot 10^{-4}$ the dropped exponential tail carries about 3% of the variance of a zero-mutation branch (0.4% of its mean).

- **Extend soft tails onto the grid** for a fixed number of decay lengths before the FFT.
- **Exponential-kernel recursion:** convolution with $p_e e^{-k t}$ against a piecewise-linear function is a first-order recursion, exact for the interpolant and $O(n)$ at any $k \cdot dx$.
- **Linearity split:** convolve a gridded core and analytic tail pieces separately and sum.

### Moment matching

Replacing an operand whose standard deviation is at most the other's spacing by a mean shift and the three-point kernel $[s/2, 1-s, s/2]$, $s = \sigma^2/dx^2$, is exact for mass, mean and variance. It applies to 0-6% of the measured convolutions, so it is not a performance lever while supports are untrimmed. Prototype with tail moments (`TailMoments`, `fn SoftTailLaw::moments()`, `fn HardApproachLaw::moments()`, `fn MassProfile::mean_variance()`, `fn MassProfile::mass_quantile()`) on the unmerged branch `moment-matched-convolution`.

## Experiment results (2026-10-19)

Prototype in the `experiment(...)` commit `6d06d211a` on branch `claude-grid-trim` (base `26c1c8533`), switched by environment variables: `TT_TRIM=c|a|p`, `TT_COARSE`, `TT_TAILS`, `TT_GRID_STATS`, `TT_DUMP_NODE`; not for merging. Trimming keeps a trimmed hard side `Hard` at the new edge (option c); option (a), a `Linear` tail fitted at the new edge, put mass past child tip dates on the simulated data (a parent's interval reached 0.15 years past its child's sampling date) and was dropped. Release build, `--jobs 1`, idle machine.

### Where to trim: the forward posterior only

| variant | what is trimmed by mass, hard sides included | result |
| --- | --- | --- |
| c | forward convolution output and forward posterior | `data/mpox/clade-ii/1000`: 152 internal nodes without a date (empty posterior), tips with uncertain dates moved by up to 6 years; `data/sc2/4500`: some intervals at 3% of base |
| p | forward posterior only | consistent with base on all datasets (below); without the spacing cap, `mpox` exceeds `--max-grid-points` |
| p + coarse | as p, convolution output spacing at most the coarser operand's | below |

Trimming the convolution output by its own mass is unsafe: the subtree distribution it is multiplied with can put the product's mass exactly where the convolution was trimmed (the multiplicative hazard of the empty invariant, [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md#empty-invariant)). The backward messages must not be trimmed for the same reason (a first prototype that trimmed every `fn convolve_across_edge()` call cut a subtree message at 2024.22 where the posterior had mass up to 2024.38). Trimming the posterior is safe in exact arithmetic: a parent's posterior mass removed below its $\epsilon$ quantile bounds the mass any child loses by $\epsilon$. Children still get small grids, because their operand, the parent's cavity, comes from the trimmed posterior.

Trimming and the spacing cap have to land together, not in the order of the [issues](#issues): a trimmed posterior is gridded at 300 points over its mass window, its fine spacing passes into the next convolution with a wide branch distribution, and the forward grids grow (sc2: median forward convolution 9 989 to 13 349 points, maximum 240 050; mpox: a 1.5-million-point grid).

### Speed (p + coarse)

| dataset | wall time base | p + coarse | forward pass, one round | peak memory |
| --- | --- | --- | --- | --- |
| validation `const/rep00/L10000_mu0.01`, fixed clock, 3 passes | 36.3 s | 5.9 s | 8.7 s to 0.3 s | 398 to 123 MB |
| `data/sc2/4500` | 47.3 s | 14.9 s | 9.8 s to 1.05 s | 1070 to 695 MB |
| `data/mpox/clade-ii/1000`, fixed clock | 5.6 s | 5.3 s | 0.33 s to 0.24 s | 414 to 410 MB |
| `data/ebola/362` | 1.0 s | 0.9 s | 0.08 s to 0.06 s | 60 MB |

Median forward grid (convolution / posterior points): simulated 85 424 / 84 245 to 459 / 461; sc2 9 989 / 9 620 to 449 / 300. On this base mpox no longer has the 30 000-point grids measured above, so it gains little. Validation suite (2880 runs): run time ratio per variant 0.50-0.72, total 1300 s to 516 s on 12 cores.

### Accuracy

- Validation suite, 2880 runs against base: no failures, no node dated before its parent, median change of date RMSE per variant below $2\cdot10^{-4}$ years. 33 runs change their date RMSE by more than 0.05 years (16 worse, 17 better), all at $L=10^3$, $\mu=3\cdot10^{-4}$, where the root date is off by 7-18 years in base too.
- Real data against base (no truth): every node dated; largest date change 0.005 years (ebola), 0.027 (mpox), 0.009 (sc2); interval widths median 0.993-0.997 of base. On sc2, 97 of 5099 internal intervals are narrower by more than 10%, all on a ladder of nodes pinned at the earliest tip date, where base itself dates a child interval 0.034 years before its parent's (one violation above 0.01 years in base, none trimmed).
- Calibration gate (240 sets, `--coalescent opt`, estimated rate, `--time-marginal only-final`): 90% coverage unchanged within 0.013 wherever $\mu \ge 10^{-3}$ or $L=10^4$. At $L=10^3$, $\mu=3\cdot10^{-4}$ coverage falls from 0.830 to 0.743, in all 30 runs, with interval widths down 20-35% and the date RMSE unchanged. At this signal base posteriors are gridded coarsely over their whole hard support (one node: 68 points over five years), and base reports a 90% interval of 1.2 years for a posterior whose 5-95% quantiles span 0.15 years; trimmed posteriors are resolved at 300 points over their mass. The lower coverage is the model's under-coverage at low signal, which the coarse base grids had masked; it needs its own fix (rate and prior at low signal), not wider grids.

### Branch distributions (issue 4)

With the forward pass fixed, building the branch distributions is half the run time of the simulated set (0.9 s of 1.9 s per round). Of that, 93% is evaluating the per-site likelihood at the 300 grid branch lengths (`fn evaluate_log_lh_on_grid()`), 5-6% building the edge contributions, about 1% re-windowing. It scales with the mutations: 0.81 s per round at $\mu=10^{-2}$, 0.11 s at $\mu=10^{-3}$ ($L=10^4$, 400 tips). Directions: fewer evaluation points with interpolation of the smooth log-likelihood, or grouping sites with the same pattern.

### Operand tails (issue 3)

Prototype: before each convolution, `Linear` tails of both operands are resampled onto the grid for 12 nats of decay at the operand's spacing (`HardApproach` tails, the branch distributions' side near zero length, left as they are). Calibration gate as above, without trimming: at $L=10^3$, $\mu=3\cdot10^{-4}$ coverage 0.830 to 0.817 and date RMSE 0.494 to 0.586; at a fixed true rate on one set the interval widths stay within 1% (median 0.994) and dates within 0.14 years, while with an estimated rate the probes of the rate step change (one set: rate $1.9\cdot10^{-4}$ to $2.3\cdot10^{-4}$, true $3\cdot10^{-4}$). The prototype also fails on 31 of 240 runs, all with $\mu \ge 10^{-3}$, $L=10^4$ or $\mu=10^{-2}$ (`NaN` when re-peaking a backward message), so its high-signal numbers are biased. No sign that the dropped tails limit the low-signal coverage; the coarse-grid effect above is the larger one. Not pursued.

### Recommendation

Land posterior trimming (option c at the posterior only) together with the spacing cap `dx = min(window / (GRID_POINTS - 1), max(dx_a, dx_b))`, both in the forward pass only, as one change; record in a decision that the backward messages and the forward convolution output keep their hard support, and why. Then take issue 4. Implemented: [kb/decisions/timetree-forward-posterior-trimming.md](../decisions/timetree-forward-posterior-trimming.md).

## Evidence from a regrid stack built before the experiments

On 2026-10-18, before the experiments above, a stack of changes on `rust` of that day replaced the grid handling of time inference with an exact mass model of piecewise log-linear distributions, an exact convolution of ranges, the v0 grid rule of the FFT convolution, and a regrid of every node distribution at mass quantiles. The stack lives on local branches only and is not merged; the commits named below are provenance, not code to read. Four of its results bear on this proposal.

### Node regrid at mass quantiles does not meet a 1% accuracy rule

- **Construction**: after the messages of a node are combined, in the backward and in the forward pass, the node distribution is resampled to $n$ points over a window. A hard side stays at its grid edge. A soft side ends at the nearer of the $\epsilon$ mass quantile and the point where the density falls to $10^{-15}$ of the peak. When the resampled grid changes the window mass by more than 1%, the node keeps its input grid
- **Study**: $\epsilon \in \{10^{-3}, 10^{-4}, 10^{-6}, 10^{-9}\}$ times $n \in \{100, 200, 300, 600, 1200\}$ on `flu/h3n2/20`, `flu/h3n2/200`, `ebola/362`, `zika/86`, `rsv/a/500`, `mpox/clade-ii/500`, `dengue/500` and `lassa/L/200`, against the same pipeline without the regrid. The accuracy rule: every node date and both 90% bounds stay within 1% of the reference 90% width
- **Result**: no setting meets the rule. The best setting, $\epsilon = 10^{-9}$ with 100 points, moves one node by 0.213 of its 90% width (`mpox/clade-ii/500`, `OL504743`); $\epsilon = 10^{-6}$ with 200 points reaches 0.532, $\epsilon = 10^{-3}$ with 300 points 3.2 (`dengue/500`). Errors do not fall as $n$ grows to 1200
- **Cause**: the window spans many 90% widths, so $n$ points leave few points inside the 90% interval, and a node date is a grid peak, so its error is of the order of the window spacing. With 5000 to 20 000 points the deviation falls to about 1% on `zika/86`, `ebola/362` and `flu/h3n2/200`
- **Fallback**: 13% to 44% of the regrid calls keep their input, fewer as $n$ grows (0% at 20 000 points)
- **Grid size**: the regrid does not bound it. The largest grid of every setting, 440 000 to 984 000 points, is a convolution grid, of the order of the reference (793 186 points on `lassa/L/200`), and some settings exceed the grid point limit of 1 000 000 where the reference does not
- **Conclusion**: a fixed point count over a mass window costs accuracy and leaves the largest grids, the convolution grids, unchanged. This proposal trims only the forward posterior and caps the spacing of the convolution output instead
- **Provenance**: local branch `worktree/feat/timetree-node-regrid-parameters`, commits `e085de2e1` and `adbbf3f65`

### Cost of the v0 convolution grid rule

v0 convolves at the spacing $\Delta t = \max(0.005\,m,\ \min(w_a, w_b) / 150)$ [packages/legacy/treetime/treetime/node_interpolator.py#L162](../../packages/legacy/treetime/treetime/node_interpolator.py#L162), where $m$ is the branch length of one mutation and $w_a$, $w_b$ are the full widths at half maximum of the two operands, and samples each operand from the start of its effective support (density above $10^{-15}$ of the peak) over twice the larger support width [packages/legacy/treetime/treetime/node_interpolator.py#L199-L209](../../packages/legacy/treetime/treetime/node_interpolator.py#L199-L209).

Release builds on the host, `-j 1`, run time by commit of the stack:

- **The v0 spacing rule** (`00a2d21da`): `ebola/362` 1.1 s to 2.1 s
- **The hard side of a cropped FFT result stored at the FFT spacing out to the operand ends** (`13fa8ddb0`): `flu/h3n2/200` 1.45 s to 3.6 s, with 3 to 6 times the memory on cases of the validation suite
- **Soft operand tails sampled past the effective support by the larger support width** (`38a325f2b`): `ebola/362` 3.2 s to 5.9 s
- **The whole stack**: 5.9 s (`ebola/362`), 10.2 s (`flu/h3n2/200`), 175 s and 2.15 GB (`sc2/4500`), against 1.0 s, 1.3 s, 33 s and 1.0 GB on `rust` of 2026-10-19. 7 of 36 runs of a validation subset failed, 6 of them on the grid point limit; on the bisected cases all failures appear with the node regrid, and without it one failure remains (a backward message gridded over \[1080.26, 2024.7\] at spacing $4.7 \cdot 10^{-4}$, 2.0 million points)

The stored extension is expensive because its length is the distance from the crop point to the physical hard end of the result, the sum of the operands' hard ends. In the forward pass that end is inherited from the root (see [Mechanism](#mechanism)), so the extension spans the whole inherited support at the FFT spacing. v0 represents the same extension with 9 points between the crop point and its bound, whatever the distance [packages/legacy/treetime/treetime/node_interpolator.py#L251-L254](../../packages/legacy/treetime/treetime/node_interpolator.py#L251-L254). A representation of constant size avoids this cost.

### A hard edge at an FFT crop point removes density that another message needs

- **Construction**: the FFT result is cropped where it falls below $10^{-13}$ of its peak, as in v0 [packages/legacy/treetime/treetime/node_interpolator.py#L226](../../packages/legacy/treetime/treetime/node_interpolator.py#L226), and the side of the result opposite its soft tail is set `Hard` at the crop point
- **Effect**: on `lassa/L/50` with `--max-iter=0`, the posterior of `NODE_0000005` peaked on its own cut edge, 1908.69 with the 90% interval \[1906.05, 1908.69\]. With the result extended by the 3-cell slope to its physical hard end, the node is dated 1853.64 \[1828.23, 1875.55\]; before the stack changed the FFT convolution (commit `58a00ff0c`) it was dated 1850.16 \[1825.48, 1871.83\]
- **v0**: v0 extends both cropped ends with the slope of the last 3 cells out to its bounds `Tmin` and `Tmax` [packages/legacy/treetime/treetime/node_interpolator.py#L241-L256](../../packages/legacy/treetime/treetime/node_interpolator.py#L241-L256)
- **Conclusion**: cutting a message where its own density is small removes density where the product with another message has its mass. This is the hazard behind the rule above that convolution outputs and backward messages are not trimmed
- **Provenance**: commit `13fa8ddb0` of the stack

### A point shortcut for a narrow branch picks the wrong mode of a two-mode node

- **Construction**: when the effective support of the branch distribution is narrower than the full width at half maximum of the node distribution divided by 150, the stack shifts the node distribution by the branch peak instead of convolving
- **Case**: a node distribution with modes at 0 and 100, standard deviations 0.001 and 0.1, heights 1 and 0.99, and a Gaussian branch distribution with standard deviation 0.02. The width at half maximum of the node is that of its wide mode, so the shortcut applies. The full convolution lowers the narrow mode to 0.04994 and the wide mode to 0.97077, so the mode moves to 100; the point shortcut keeps both heights and the mode at 0
- **Conclusion**: a gate on the global width of the node has no error bound, because a node feature narrower than the branch passes it. At a node feature with standard deviation $\sigma_f$, the shortcut errs by about $\sigma_b^2 / (2\sigma_f^2)$ in neg-log when $\sigma_b \ll \sigma_f$, where $\sigma_b$ is the standard deviation of the branch distribution, and by more as $\sigma_b$ approaches $\sigma_f$. The [moment matching](#moment-matching) of this proposal gates on the other operand's grid spacing instead: a grid that resolves a feature narrower than the branch has a spacing below $\sigma_b$, so the branch is convolved in full
- **Provenance**: test `test_convolution_narrow_branch_point_misplaces_the_mode_of_a_two_mode_node` in `packages/treetime-distribution/src/distribution_ops/__tests__/test_convolve_function.rs` at commit `38a325f2b` (lines 190 to 245); the file does not exist on `rust`

## Issues

1. Trim hard sides by mass: implemented ([kb/decisions/timetree-forward-posterior-trimming.md](../decisions/timetree-forward-posterior-trimming.md))
2. Output spacing: implemented with 1
3. [M-distribution-convolution-drops-operand-tails.md](../issues/M-distribution-convolution-drops-operand-tails.md): operand tails in the convolution
4. [N-timetree-branch-grid-evaluated-one-length-at-a-time.md](../issues/N-timetree-branch-grid-evaluated-one-length-at-a-time.md): the next bottleneck after 1 and 2

## Related

- [M-timetree-node-dates-move-when-a-time-grid-gains-one-point.md](../issues/M-timetree-node-dates-move-when-a-time-grid-gains-one-point.md): node dates at the grid argmax move when a grid changes by one point
- [M-timetree-near-flat-soft-tail-stretches-the-time-grid.md](../issues/M-timetree-near-flat-soft-tail-stretches-the-time-grid.md): a soft tail with a slope near zero puts a mass quantile far outside the grid
- [M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md](../issues/M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md): a narrow branch distribution passes its spacing to the forward product, which then needs millions of points
- [M-distribution-fft-convolution-hard-edge-first-order-error.md](../issues/M-distribution-fft-convolution-hard-edge-first-order-error.md): first-order error of the FFT convolution at a hard edge
- [exponential-branch-length-tails.md](exponential-branch-length-tails.md)
