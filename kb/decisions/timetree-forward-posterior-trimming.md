# Forward-pass posteriors are trimmed to their probability mass

**Type**: Intentional divergence from v0, and a numerical approximation of the marginal time inference.

v1 cuts every refined node-time posterior of the forward pass at its lower and upper $\epsilon_\text{trim}$ mass quantiles, with $\epsilon_\text{trim} = 10^{-9}$ (`TRIM_EPS`), and caps the grid spacing of forward convolutions. Backward messages and forward convolution outputs keep their full support. v0 has no mass trim: it grids each distribution by its own interpolation rules ([packages/legacy/treetime/treetime/node_interpolator.py#L162-L256](../../packages/legacy/treetime/treetime/node_interpolator.py#L162-L256)).

## Problem

Before the trim, every forward posterior kept a hard lower edge inherited from the root, and every forward convolution output kept the finest spacing of any factor upstream. Forward grids held a median 17 to 124 times more points than the interval that holds their mass, all on the hard left side, where the density at the grid edge is below $10^{-50}$ of the peak. Forward convolution, division and multiplication took about 95% of the distribution-operation time. Measurements and mechanism: [kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md).

## Rule

### What is trimmed

`fn refine_distribution_from_parent()` [packages/treetime/src/timetree/inference/forward_pass.rs#L170-L221](../../packages/treetime/src/timetree/inference/forward_pass.rs#L170-L221) forms the refined posterior of a node, the normalized product of the forward message and the node's subtree distribution, and trims it with `fn ResolvedForward::trim_posterior_to_mass()` [packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L151-L168](../../packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L151-L168).

- **Trimmed**: the refined posterior of each internal node and of each leaf without an exact date
- **Not trimmed**: the root posterior, exactly dated nodes, a node without a subtree distribution (its posterior is the forward message), every backward message, and every forward convolution output
- **Left unchanged**: `Empty`, `Point` and `Range` posteriors, a function without a finite peak, and a function whose total mass the mass model cannot size (an `Error` side, or a soft tail that does not decay). The contradiction check after the trim then still sees an empty product and keeps the node's given date

### Where the cuts lie

`fn trim_with_profile()` [packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L54-L91](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L54-L91) computes the quantiles on the mass model of the untrimmed posterior ([timetree-posterior-statistics-count-tail-mass.md](timetree-posterior-statistics-count-tail-mass.md)):

- **Grid**: the lower cut snaps down to the largest grid point at or below the lower $\epsilon_\text{trim}$ quantile, the upper cut up to the smallest grid point at or above the upper one, and the grid arrays are sliced between them. The trim keeps the grid points and values of the input and does not resample, so the trapezoid quadrature of the samples stays unchanged. At most $\epsilon_\text{trim}$ of the mass is removed on each side
- **Hard cells**: a cut never splits the fit points of a retained hard cell (`fn GridFn::hard_cell_fit_range()`), so the power law of that cell evaluates as before
- **Cut hard side**: a `Hard` or `HardApproach` side that is cut becomes `Hard` at the cut
- **Quantile inside a `HardApproach` law**: when the quantile lies beyond the stored grid, inside the law between the grid edge and the hard bound, that side is not cut and keeps its law. A slice cannot end inside a tail law, and a cut at the grid edge would drop the law's mass (for density 1 on $[1, 2]$ with a flat law down to 0, half the mass)
- **Cut soft side**: a `Linear` side keeps a `Linear` law that holds exactly the mass the cut removed. With $p(t_c)$ the density at the cut point $t_c$ and $m_\text{beyond}$ the model mass beyond it, both relative to the peak, the new neg-log slope has the magnitude $p(t_c) / m_\text{beyond}$ (`fn tail_with_mass_beyond()` [packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L252-L267](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L252-L267))

### Spacing cap of forward convolutions

The forward convolution resamples its output onto its mass window at the spacing

$$\Delta t = \min\left(\frac{w}{N - 1},\ \max(\Delta t_a, \Delta t_b)\right)$$

where:

- $w$ -- width of the mass window of the convolution output
- $N$ -- `GRID_POINTS`, 300
- $\Delta t_a$, $\Delta t_b$ -- grid spacings of the two function operands

The backward convolution keeps the finer of the window spacing and the FFT spacing. Without the cap, a trimmed posterior with 300 points over its mass window passes its fine spacing into the next convolution with a wide branch distribution, and the forward grids grow (sc2: median forward convolution output 9 989 to 13 349 points, largest 240 050; mpox: a grid of 1.5 million points).

A capped message can be too coarse for the posterior it forms when the subtree distribution is much narrower than the message. `fn ForwardConvolution::resolve_posterior()` [packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L51-L141](../../packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L51-L141) therefore checks each posterior: when it has fewer than $N - 1$ cells across its $\epsilon_\text{trim}$ mass window and a finer spacing exists, the message is resampled once at $\max(\Delta t_\text{fine}, w_p / (N - 1))$, where $\Delta t_\text{fine}$ is the FFT spacing and $w_p$ the width of the posterior's mass window; when that posterior is still under-resolved, the message takes the FFT spacing. The trim reuses the mass model of this check.

This check replaced an acceptance test on one run that cannot decide it: on `data/flu/h3n2/200` at `--clock-rate=3.058e-3 --keep-root`, capping a single node without a subtree distribution moves `NODE_0000138` by 0.35 of its 90% width, and the uncapped 300-point run lies 0.55 to 1.4 widths away from the converged references on the cluster nodes ([kb/issues/H-timetree-backward-messages-too-heavy-far-below-their-peak.md](../issues/H-timetree-backward-messages-too-heavy-far-below-their-peak.md)).

### Confidence intervals are clipped at physical bounds

`NodePosterior.subtree_support` stores the hard sides of the node's subtree distribution (`fn Distribution::support_bounds()`), and `fn node_confidence_interval()` [packages/treetime/src/timetree/confidence.rs#L113-L135](../../packages/treetime/src/timetree/confidence.rs#L113-L135) clips the interval to them. Neither the node's own cut nor an ancestor's cut ends an interval: with `--confidence`, the rate-uncertainty term can widen an interval past the trimmed support. The posterior before the trim is not a physical support either, because the cavity keeps the hard edge of its dividend, the trimmed parent posterior.

## Why trimming is safe, and where it is not

The joint density of a parent time and a child time integrates over the child time to the parent's posterior. Removing the parent's posterior mass beyond its $\epsilon_\text{trim}$ quantiles therefore removes at most $\epsilon_\text{trim}$ per side of every child's mass, and along a root-to-tip path at most $2\,\epsilon_\text{trim}$ per trimmed ancestor. A child's posterior cannot become empty from the trim, because the parent's posterior holds the child's subtree evidence.

This bound assumes that the backward messages and the forward joint agree. Far below a message peak they do not: the FFT result is trusted only above $10^{-13}$ of its peak, and below that the message is an extrapolated tail law that can be hundreds of nats too heavy. On `data/flu/h3n2/200` at a fixed rate a child then holds 23% of its mass below its parent's cut, where the parent holds $6.4 \cdot 10^{-6}$, which exact arithmetic forbids. A cut of the parent at its quantile then removes mass that the child needs. On the same data, the trim moves `NODE_0000064`, `065` and `084` by about 1 year toward earlier dates, away from the 9600-point reference; with the trim off they stay within their 90% width. The untrimmed dates of that region are not reliable either: references at 4800 and 9600 points differ by 0.6 to 3 years on these nodes ([kb/issues/H-timetree-backward-messages-too-heavy-far-below-their-peak.md](../issues/H-timetree-backward-messages-too-heavy-far-below-their-peak.md)).

## Why backward messages and convolution outputs keep their support

- **Convolution outputs**: the subtree distribution that the output is multiplied with can put the product's mass where the output was cut. Trimming the forward convolution output by its own mass left 152 internal nodes of `data/mpox/clade-ii/1000` without a date (empty posteriors), and moved tips with uncertain dates by up to 6 years
- **Backward messages**: a prototype that trimmed every edge convolution cut a subtree message at 2024.22 where the posterior had mass up to 2024.38. A hard edge set where a message's own density is small removes density where the product with another message has its mass (proposal, section "A hard edge at an FFT crop point removes density that another message needs")

## Rejected alternatives

- **`Linear` tail fitted at the cut**: put mass past child tip dates (a parent's interval reached 0.15 years past its child's sampling date)
- **`Linear` tail truncated at the original hard bound**: needs a new `BoundaryBehavior` variant in multiplication, division, mass and tail composition, and the next forward convolution, which reads grid values only, loses the gain ([kb/issues/M-distribution-convolution-drops-operand-tails.md](../issues/M-distribution-convolution-drops-operand-tails.md))
- **Resample the trimmed posterior onto 300 points**: on a coarse posterior the log-linear interpolation adds no information but changes the quadrature from the trapezoid of the samples to the integral of the log-linear interpolant. In the Kingman test against the exact posterior (`packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs`) it doubled the largest CDF error (0.0065 to 0.0114) at every tested $\epsilon_\text{trim}$, while the spacing cap alone left it unchanged
- **The v0 convolution spacing** $\max(0.005\,m,\ \min(w_a, w_b) / 150)$, where $m$ is the branch length of one mutation and $w_a$, $w_b$ the widths at half maximum of the operands: 2 to 3 times the run time (proposal, section "Cost of the v0 convolution grid rule")

## Choice of $\epsilon_\text{trim}$

### Selection rule

The value is the smallest of $\{5 \cdot 10^{-4}, 10^{-5}, 10^{-7}, 10^{-9}\}$ that meets the rule below; the run time is flat across these values. The untrimmed run is not an oracle: at low signal its posteriors have 68 points over five years. Each run is compared with a converged reference: the untrimmed pipeline without the cap at 4800 points, whose nodes count as converged when a further doubling to 9600 points moves their date and 90% bounds by at most 1% of their 90% width. For each converged node with reference 90% width $W > 0$:

$$e = \frac{\max(|\Delta d|, |\Delta l|, |\Delta u|)}{W}$$

where:

- $\Delta d$, $\Delta l$, $\Delta u$ -- differences of the node date, the lower and the upper 90% bound to the reference

A value passes a run when (1) the run succeeds, dates every node that the untrimmed 300-point run dates, and violates no parent-child edge that the untrimmed run does not violate, and (2) the 99th percentile of $e$ is at most that of the untrimmed run plus 0.01 and the largest $e$ at most the untrimmed largest plus 0.05.

Runs: `data/ebola/362`, `data/mpox/clade-ii/500`, `data/sc2/5100` and `data/flu/h3n2/200` at the fixed rate that a default run estimates, with `--keep-root --time-marginal always --seed=1`; the simulated sets `const/rep00/L1000_mu0.0003` at its true rate (val-lo) and with `--coalescent opt` and an estimated rate (val-lo-est). val-lo-est is judged by its fixed-rate twin, because its rate changes with the setting.

### Outcome at $\epsilon_\text{trim} = 10^{-9}$

| Run | Converged nodes | p99 / largest $e$, untrimmed | p99 / largest $e$, $10^{-9}$ | Result |
| --- | --- | --- | --- | --- |
| ebola/362 | 343 of 359 | 0.0543 / 0.0594 | 0.0543 / 0.0596 | pass |
| val-lo | 399 of 399 | 0.0890 / 0.0913 | 0.0869 / 0.0896 | pass |
| val-lo-est | 287 of 399 | 0.5208 / 0.5979 | 0.3989 / 0.5037 | pass |
| mpox/clade-ii/500 | 594 of 618 | 0.0209 / 0.0383 | 0.0209 / 0.0435 | pass |
| flu/h3n2/200 | 32 of 196 | 0.4408 / 0.5207 | 0.3807 / 0.4402 | pass |
| sc2/5100, reference at 2400 points | 4825 of 5099 | 0.1009 / 0.3155 | 0.1101 / 0.3412 | pass |
| sc2/5100, reference at 4800 points | 4825 of 5099 | 0.1032 / 0.3210 | 0.1152 / 0.3448 | fails p99 by 0.0020 |

Every run meets rule (1): it succeeds, dates every node that the untrimmed run dates, and violates no new parent-child edge. The maintainers accepted $10^{-9}$ with the sc2 residual on 2026-10-20: the p99 of $e$ exceeds its limit by 0.0020 against the 4800-point reference and passes against the 2400-point reference, while the two references differ from each other by a p99 of 0.032. No sc2 reference at 9600 points exists, so the converged nodes of sc2 are those that agree between 2400 and 4800 points. $10^{-7}$ passes sc2 against both references and fails flu by 0.008.

Flu passes because the posterior check of the spacing cap gives under-resolved posteriors a finer message (see above); its largest forward grid is 67 449 points (FFT of 147 456) against 7 876 points untrimmed, and the run still takes 1.2 s of CPU. With the cap off and the trim kept, flu fails (p99 0.504); with the trim off and the cap kept, flu (p99 0.554) and val-lo-est fail.

Cost on one host, CPU time and peak memory: `data/sc2/5100` untrimmed 73.3 s and 953 MB, at $10^{-9}$ 21.7 s and 885 MB, at $10^{-7}$ 21.4 s and 868 MB; with `--coalescent opt` 78.3 s and 1004 MB untrimmed, 26.2 s and 995 MB at $10^{-9}$; the simulated set `const/rep00/L10000_mu0.01` 56.2 s and 188 MB untrimmed, 9.0 s and 171 MB at $10^{-9}$.

## Results

The measurements of this file come from outputs equal to those of commit `6868c7087`. A later correction of the mass model (each cell takes the width between its rounded end coordinates) moves dates of the measured ebola, flu and zika runs by at most $2.3 \cdot 10^{-10}$ of a 90% width, and on `data/sc2/5100` by up to 0.042 of a width: the first difference is a one-ulp move of the parent of `USA/MN-MDH-48906/2025`, and the later rounds spread it ([kb/issues/M-timetree-node-dates-move-when-a-time-grid-gains-one-point.md](../issues/M-timetree-node-dates-move-when-a-time-grid-gains-one-point.md)).

Release builds on one host, one thread, CPU time (user plus system) and peak resident memory, medians of three alternating rounds, against a release build of `rust` `d86d3eb65` (before this change and the mass model of [timetree-posterior-statistics-count-tail-mass.md](timetree-posterior-statistics-count-tail-mass.md)), measured on the final code including the later mass-model correction:

| Run | CPU time | CPU ratio | Peak memory | Memory ratio |
| --- | --- | --- | --- | --- |
| `data/sc2/5100` | 48.64 to 19.30 s | 0.397 | 665 to 622 MB | 0.936 |
| `data/ebola/362` | 1.17 to 1.13 s | 0.966 | 56 to 52 MB | 0.927 |
| `data/flu/h3n2/200` | 1.94 to 0.70 s | 0.361 | 35 to 33 MB | 0.940 |
| 36-run validation subset | 351.1 to 122.5 s (total) | 0.349 | 123 to 133 MB (largest run) | 1.084 |

No validation run takes more than 1.11 times its earlier CPU time (1.39 to 1.54 s).

**Accepted residual in peak memory**: 12 of the 36 validation runs have more than 1.1 times their earlier peak memory, at most 1.193 (`keep-opt-est` rep 1, $L = 10^3$, $\mu = 3 \cdot 10^{-4}$: 42 to 50 MB). All are small coalescent runs of 42 to 60 MB that grow by 4 to 10 MB; runs near 1.10 cross that limit between measurements. The live heap at the peak stays within about 3 MB of the earlier code. The cause is the page reuse of the jemalloc allocator, which maps new pages instead of reusing freed ones in the marginal clock-rate step; an allocator setting removes it ([kb/issues/N-alloc-jemalloc-keeps-freed-pages-per-arena.md](../issues/N-alloc-jemalloc-keeps-freed-pages-per-arena.md)). The maintainers accepted this residual on 2026-10-20 and keep the allocator setting in that issue.

Accuracy on the validation subset (36 runs with known truth, `keep-none-fixed`, `keep-opt-est`, `reroot-opt-est` at $L \in \{10^3, 10^4\}$, $\mu \in \{3 \cdot 10^{-4}, 10^{-3}, 10^{-2}\}$): no failed run, all 28 764 nodes dated; 16 of 18 cells within $\max(1.1 \cdot \text{base}, \text{base} + 0.05)$ years of root error. The two cells above the limit (`keep-opt-est` at $L = 10^4$, $\mu = 3 \cdot 10^{-4}$: 0.528 to 0.585 years, limit 0.581; `reroot-opt-est` at $L = 10^3$, $\mu = 10^{-3}$: 1.828 to 2.138 years, limit 2.011) come from the clock rate step on 300-point messages, not from the trim ([kb/issues/H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md](../issues/H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md)).

Calibration (12 sets, `keep-opt-est` with `--time-marginal only-final`, mean of 2 replicates): 90% coverage of node dates falls by at most 0.024 in any cell, at $L = 10^3$, $\mu = 3 \cdot 10^{-4}$ (0.807 to 0.783, mean interval width 0.400 to 0.353 years). `--confidence` only widens intervals: per cell it raises the coverage by 0 to 0.0038 and the mean width by 0.0012 to 0.0024 years, and no single run changes by more than 0.0050 in coverage and 0.0033 years in width.

| Cell | 90% coverage | Mean width (years) | 90% coverage with `--confidence` | Mean width with `--confidence` |
| --- | --- | --- | --- | --- |
| $L = 10^3$, $\mu = 3 \cdot 10^{-4}$ | 0.783 | 0.353 | 0.787 | 0.355 |
| $L = 10^3$, $\mu = 10^{-3}$ | 0.867 | 0.376 | 0.867 | 0.377 |
| $L = 10^3$, $\mu = 10^{-2}$ | 0.899 | 0.243 | 0.899 | 0.245 |
| $L = 10^4$, $\mu = 3 \cdot 10^{-4}$ | 0.848 | 0.316 | 0.850 | 0.318 |
| $L = 10^4$, $\mu = 10^{-3}$ | 0.866 | 0.228 | 0.867 | 0.231 |
| $L = 10^4$, $\mu = 10^{-2}$ | 0.905 | 0.093 | 0.907 | 0.094 | An earlier prototype that resampled the trimmed posteriors lowered the coverage at the lowest signal from 0.830 to 0.743; slicing without resampling does not.

## Code references

- `const TRIM_EPS` [packages/treetime/src/timetree/inference/runner.rs#L40](../../packages/treetime/src/timetree/inference/runner.rs#L40)
- `fn refine_distribution_from_parent()` [packages/treetime/src/timetree/inference/forward_pass.rs#L170-L221](../../packages/treetime/src/timetree/inference/forward_pass.rs#L170-L221)
- `fn convolve_forward_across_edge()`, `struct ForwardConvolution`, `fn ResolvedForward::trim_posterior_to_mass()` and `fn ResolvedForward::into_untrimmed()`, both returning `struct ForwardPosterior` [packages/treetime-distribution/src/distribution_ops/edge_convolution.rs](../../packages/treetime-distribution/src/distribution_ops/edge_convolution.rs)
- `fn trim_to_mass()`, `fn trim_with_profile()`, `enum SpacingLimit` [packages/treetime-distribution/src/distribution_ops/mass_domain.rs](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs)
- `NodePosterior.subtree_support` [packages/treetime/src/timetree/inference/result.rs#L103](../../packages/treetime/src/timetree/inference/result.rs#L103)

## Related

- [kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md): measurements, prototype and recommendation
- [distribution-tails-and-arithmetic.md](distribution-tails-and-arithmetic.md): tail assignments of the timetree messages
- [timetree-posterior-statistics-count-tail-mass.md](timetree-posterior-statistics-count-tail-mass.md): the mass model that places the cuts
