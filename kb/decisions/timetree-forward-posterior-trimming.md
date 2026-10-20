# Forward-pass posteriors are trimmed to their probability mass

**Type**: Intentional divergence from v0, and a numerical approximation of the marginal time inference.

v1 cuts every refined node-time posterior of the forward pass at its lower and upper $\epsilon_\text{trim}$ mass quantiles, with $\epsilon_\text{trim} = 10^{-9}$ (`TRIM_EPS`), and caps the grid spacing of forward convolutions. Backward messages and forward convolution outputs keep their full support. v0 has no mass trim: it grids each distribution by its own interpolation rules ([packages/legacy/treetime/treetime/node_interpolator.py#L162-L256](../../packages/legacy/treetime/treetime/node_interpolator.py#L162-L256)).

## Problem

Before the trim, every forward posterior kept a hard lower edge inherited from the root, and every forward convolution output kept the finest spacing of any factor upstream. Forward grids held a median 17 to 124 times more points than the interval that holds their mass, all on the hard left side, where the density at the grid edge is below $10^{-50}$ of the peak. Forward convolution, division and multiplication took about 95% of the distribution-operation time. Measurements and mechanism: [kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md).

## Rule

### What is trimmed

`fn refine_distribution_from_parent()` [packages/treetime/src/timetree/inference/forward_pass.rs#L170-L222](../../packages/treetime/src/timetree/inference/forward_pass.rs#L170-L222) forms the refined posterior of a node, the normalized product of the forward message and the node's subtree distribution, and trims it with `fn ResolvedForward::trim_posterior_to_mass()` [packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L155-L174](../../packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L155-L174).

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

A capped message can be too coarse for the posterior it forms when the subtree distribution is much narrower than the message. `fn ForwardConvolution::resolve_posterior()` [packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L51-L140](../../packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L51-L140) therefore checks each posterior: when it has fewer than $N - 1$ cells across its $\epsilon_\text{trim}$ mass window and a finer spacing exists, the message is resampled once at $\max(\Delta t_\text{fine}, w_p / (N - 1))$, where $\Delta t_\text{fine}$ is the FFT spacing and $w_p$ the width of the posterior's mass window; when that posterior is still under-resolved, the message takes the FFT spacing. The trim reuses the mass model of this check.

This check replaced an acceptance test on one run that cannot decide it: on `data/flu/h3n2/200` at `--clock-rate=3.058e-3 --keep-root`, capping a single node without a subtree distribution moves `NODE_0000138` by 0.35 of its 90% width, and the uncapped 300-point run lies 0.55 to 1.4 widths away from the converged references on the cluster nodes ([kb/issues/H-timetree-backward-messages-too-heavy-far-below-their-peak.md](../issues/H-timetree-backward-messages-too-heavy-far-below-their-peak.md)).

### Confidence intervals are clipped at physical bounds

`NodePosterior.subtree_support` stores the hard sides of the node's subtree distribution (`fn Distribution::support_bounds()`), and `fn node_confidence_interval()` [packages/treetime/src/timetree/confidence.rs#L113-L135](../../packages/treetime/src/timetree/confidence.rs#L113-L135) clips the interval to them. Neither the node's own cut nor an ancestor's cut ends an interval: with `--confidence`, the rate-uncertainty term can widen an interval past the trimmed support. The posterior before the trim is not a physical support either, because the cavity keeps the hard edge of its dividend, the trimmed parent posterior.

## Why trimming is safe, and where it is not

The joint density of a parent time and a child time integrates over the child time to the parent's posterior. Removing the parent's posterior mass beyond its $\epsilon_\text{trim}$ quantiles therefore removes at most $\epsilon_\text{trim}$ per side of every child's mass, and along a root-to-tip path at most $2\,\epsilon_\text{trim}$ per trimmed ancestor. A child's posterior cannot become empty from the trim, because the parent's posterior holds the child's subtree evidence.

This bound assumes that the backward messages and the forward joint agree. Far below a message peak they do not: the FFT keeps about 37 nats of dynamic range, and below that the message is an extrapolated tail law that can be hundreds of nats too heavy. On `data/flu/h3n2/200` at a fixed rate a child then holds 23% of its mass below its parent's cut, where the parent holds $6.4 \cdot 10^{-6}$, which exact arithmetic forbids. A cut of the parent at its quantile then removes mass that the child needs. The untrimmed dates of that region are not reliable either: references at 4800 and 9600 points differ by 3 to 8 years there ([kb/issues/H-timetree-backward-messages-too-heavy-far-below-their-peak.md](../issues/H-timetree-backward-messages-too-heavy-far-below-their-peak.md)).

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

Runs: `data/ebola/362`, `data/mpox/clade-ii/500`, `data/sc2/4500` and `data/flu/h3n2/200` at the fixed rate that a default run estimates, with `--keep-root --time-marginal always --seed=1`; the simulated sets `const/rep00/L1000_mu0.0003` at its true rate (val-lo) and with `--coalescent-opt` and an estimated rate (val-lo-est). val-lo-est is judged by its fixed-rate twin, because its rate changes with the setting.

### Outcome at $\epsilon_\text{trim} = 10^{-9}$

| Run | Converged nodes | p99 / largest $e$, untrimmed | p99 / largest $e$, $10^{-9}$ | Result |
| --- | --- | --- | --- | --- |
| ebola/362 | 343 of 359 | 0.0543 / 0.0594 | 0.0465 / 0.0512 | pass |
| val-lo | 399 of 399 | 0.0890 / 0.0913 | 0.0891 / 0.0928 | pass |
| val-lo-est | 312 of 399 | 0.5366 / 0.6828 | 0.5267 / 0.6684 | pass |
| mpox/clade-ii/500 | 594 of 618 | 0.0209 / 0.0383 | 0.0226 / 0.0415 | pass |
| sc2/4500, reference at 4800 points | 4825 of 5099 | 0.1032 / 0.3210 | 0.1145 / 0.3579 | fails p99 by 0.0013 |
| flu/h3n2/200, cap off | 32 of 196 | 0.4408 / 0.5207 | 0.504 (p99) | fails p99 by 0.053 |

The maintainers accepted $10^{-9}$ with these two residuals on 2026-10-20:

- **sc2**: p99 of $e$ rises by 0.011 of the 90% width against the limit 0.01. The 2400-point reference itself differs from the 4800-point reference by a p99 of 0.032, and the 4800-point reference costs 1851 s of CPU and 7.9 GB
- **flu at a fixed rate**: p99 of $e$ rises by 0.063 against the limit 0.01, on 32 converged nodes of 196. The cause is the deep-tail error of the backward messages described above, which moves the untrimmed result as well. With default flags no node of `data/flu/h3n2/200` moves by more than 0.02 of its 90% width against `rust` before the trim

The cost of the trim does not depend on $\epsilon_\text{trim}$ within the measurement noise: on `data/sc2/4500` the untrimmed run takes 164 s of CPU and 1022 MB, the run at $10^{-9}$ 41 s and 664 MB (on a loaded host).

## Results

Release builds, one thread, CPU time (user plus system) and peak resident memory, median of three alternating pairs, as ratios to the code before this change and the mass model of [timetree-posterior-statistics-count-tail-mass.md](timetree-posterior-statistics-count-tail-mass.md):

| Run | CPU time | Peak memory |
| --- | --- | --- |
| `data/sc2/4500` | 0.409 | 1.060 |
| `data/ebola/362` | 1.017 | 1.045 |
| `data/flu/h3n2/200` | 0.385 | 0.994 |

On the 36 runs of the validation subset the total CPU time fell from 351.9 s to 109.9 s (0.312) and the largest peak memory from 414 MB to 126 MB, measured on an earlier state of the same change.

Accuracy on the validation subset (36 runs with known truth, `keep-none-fixed`, `keep-opt-est`, `reroot-opt-est` at $L \in \{10^3, 10^4\}$, $\mu \in \{3 \cdot 10^{-4}, 10^{-3}, 10^{-2}\}$): no failed run, all 28 764 nodes dated; 16 of 18 cells within $\max(1.1 \cdot \text{base}, \text{base} + 0.05)$ years of root error. The two cells above the limit (`keep-opt-est` at $L = 10^4$, $\mu = 3 \cdot 10^{-4}$: 0.528 to 0.615 years; `reroot-opt-est` at $L = 10^3$, $\mu = 10^{-3}$: 1.828 to 2.043 years) come from the clock rate step, not from the trim ([kb/issues/H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md](../issues/H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md)).

Calibration (12 sets, `keep-opt-est` with `--time-marginal only-final`, mean of 2 replicates): 90% coverage of node dates falls by at most 0.010 ($L = 10^3$, $\mu = 3 \cdot 10^{-4}$: 0.807 to 0.797, mean interval width 0.400 to 0.363 years), and by at most 0.003 where $\mu \ge 10^{-3}$ or $L = 10^4$. With `--confidence` the coverage changes by at most 0.004 and the width by at most 0.003. An earlier prototype that resampled the trimmed posteriors lowered the coverage at the lowest signal from 0.830 to 0.743; slicing without resampling does not.

## Code references

- `const TRIM_EPS` [packages/treetime/src/timetree/inference/runner.rs#L40](../../packages/treetime/src/timetree/inference/runner.rs#L40)
- `fn refine_distribution_from_parent()` [packages/treetime/src/timetree/inference/forward_pass.rs#L170-L222](../../packages/treetime/src/timetree/inference/forward_pass.rs#L170-L222)
- `fn convolve_forward_across_edge()`, `struct ForwardConvolution`, `fn ResolvedForward::trim_posterior_to_mass()` [packages/treetime-distribution/src/distribution_ops/edge_convolution.rs](../../packages/treetime-distribution/src/distribution_ops/edge_convolution.rs)
- `fn trim_to_mass()`, `fn trim_with_profile()`, `enum SpacingLimit` [packages/treetime-distribution/src/distribution_ops/mass_domain.rs](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs)
- `NodePosterior.subtree_support` [packages/treetime/src/timetree/inference/result.rs#L103](../../packages/treetime/src/timetree/inference/result.rs#L103)

## Related

- [kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md): measurements, prototype and recommendation
- [distribution-tails-and-arithmetic.md](distribution-tails-and-arithmetic.md): tail assignments of the timetree messages
- [timetree-posterior-statistics-count-tail-mass.md](timetree-posterior-statistics-count-tail-mass.md): the mass model that places the cuts
