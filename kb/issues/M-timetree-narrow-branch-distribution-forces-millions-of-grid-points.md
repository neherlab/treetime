# A convolution of a narrow and a wide time distribution needs millions of grid points

When one operand of a time convolution is far narrower than the other, the convolution and the products after it grid the whole support of the wide operand at the spacing of the narrow one, and the run stops at the grid point limit ([kb/decisions/distribution-grid-point-limit.md](../decisions/distribution-grid-point-limit.md)). Two cases are known: a narrow branch-time distribution next to a wide node distribution, and a narrow node posterior next to a long branch.

## Narrow branch next to a wide node distribution

### Reproduction

```bash
treetime timetree --tree=data/mpox/clade-ii/20/tree.nwk --metadata=data/mpox/clade-ii/20/metadata.tsv \
  --metadata-id-columns=accession --alignment=data/mpox/clade-ii/20/aln.fasta.xz \
  --clock-rate=0.5 --max-iter=1 --output-all=tmp/grid-point-limit
```

The smoke matrix declares this run as the expected failure `variant = "grid-point-limit"` in `dev/smoke.toml`.

- **`rust` of 2026-10-19**: the run stops with "A grid over [2019.14, 2022] with spacing 2.17304e-7 needs 13175607 points, more than the limit of 1000000". The backtrace shows the forward pass: `fn refine_distribution_from_parent()` calls `fn distribution_multiplication()`, the product of the forward message and the subtree distribution of the first refined node
- **The forward-posterior trimming prototype** (commit `6d06d211a` on branch `claude-grid-trim`, `TT_TRIM=p TT_COARSE=1`, see [kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md)): the run stops in the same product with 13 869 392 points over \[2019.14, 2022.45\], and before that warns that `PP852952` has an empty time distribution and gets no date. Trimming the posterior does not help, because the large grid is a product of untrimmed messages
- **With the forward-posterior trim and the spacing cap** ([kb/decisions/timetree-forward-posterior-trimming.md](../decisions/timetree-forward-posterior-trimming.md), code of 2026-10-20): `data/mpox/clade-ii/1000` with `--metadata-id-columns=accession --clock-rate=0.5 --max-iter=1` stops with a grid over \[2019.48, 2023.14\] at spacing $2.04 \cdot 10^{-7}$ years, 17 968 319 points. Which operation requests this grid is not recorded
- **Default `timetree` run on `data/lassa/L/20`**: `--branch-length-mode auto` resolves to `input` (longest branch 0.49), and input mode keeps the weak positive regression rate $4.97 \cdot 10^{-5}$ ($r^2 = 0$). The run stops with "grid of 1208376 points" at the spacing 0.0102 years of a month-precision tip date, over a support of more than 11 000 years (on `rust` of 2026-10-24: 1 052 080 points over \[-8802.91, 1976.25\]). The data carry no temporal signal: the root-to-tip correlation of `treetime clock` is $r = 0.023$. The error message names weak temporal signal and `--clock-rate` as the remedy, so this is the intended guard on such input, not a wrong result. Marginal mode with `--clock-rate=4.971e-5` stops the same way; the default marginal run starts from the Gaussian rate $1.74 \cdot 10^{-4}$ and completes. v0 completes in input mode on its non-uniform grid of 125 points, with a root date of -9107

### Mechanism

- **Spacing source**: at `--clock-rate=0.5`, about 5000 times a realistic mpox rate, a branch-time distribution spans about $6.5 \cdot 10^{-5}$ years on 300 points, a spacing of $2.2 \cdot 10^{-7}$ years
- **Convolution**: `fn convolution_function_function_fine()` [packages/treetime-distribution/src/distribution_ops/convolve.rs#L271-L281](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L271-L281) resamples both operands to the finer spacing, and `fn resample_to_mass_window()` [packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L332-L368](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L332-L368) keeps that spacing for a backward result; a forward result keeps it when its posterior needs it. Convolving a node distribution several years wide at this spacing needs more points than the limit (15 583 861 on an earlier base, recorded in the decision file); when the operand is narrow, the forward message is narrow but keeps the fine spacing
- **Product**: multiplication grids the support intersection of its operands at the finest operand spacing (grid intersection contract, [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md#grid-intersection-contract)). The forward message has a soft right tail and the subtree distribution a soft left tail, so the intersection reaches from the hard left edge of the forward message to the hard right edge of the subtree distribution, 2.9 years here, at the spacing of the branch
- **Regrid stack**: on a stack of 2026-10-18 that replaced the convolution spacing with the v0 rule (`max(0.005 * one_mutation, min(FWHM_a, FWHM_b) / 150)`), the convolution needed at most 11 016 points, and the forward product of two messages 3.3 years apart needed 22 114 172 points
- **v0**: with `--time-marginal always` v0 fails on this input with an empty product ("Unexpected behavior detected in multiply function when determining peak of function with y-values '[]'", recorded in [kb/decisions/distribution-grid-point-limit.md](../decisions/distribution-grid-point-limit.md))

## Narrow node posterior next to a long branch

Four Nextstrain production builds stop in the forward pass at the fixed clock rates of their workflows, which v0 completes:

| Build | Fixed clock rate | Grid points | Grid span (years) | Stage |
| --- | ---: | ---: | --- | --- |
| rubella E1 | 0.0006 | 1 011 265 | 0 to 115.5 | time tree |
| oropouche S | 0.0014 | 1 016 064 | 0 to 7.93 | time tree |
| H5N1 cattle outbreak, PB2 | 0.00287 | 1 004 720 | 0 to 1.15 | time tree |
| H5N1 cattle outbreak, HA | 0.00455 | 1 018 346 | 0 to 1.81 | time tree at the lower rate of the rate susceptibility analysis |

- **Request**: `fn WindowedConvolution::new()` calls `fn distribution_convolution_fine()` [packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L229](../../packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L229), which resamples the branch distribution to the spacing of the parent posterior through `fn resample_dx()` [packages/treetime-distribution/src/distribution_ops/convolve.rs#L332](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L332). The span of the grid comes from the branch (7.93 years for a branch without substitutions in oropouche S, 115.5 years in rubella E1), its spacing from the posterior ($7.8 \cdot 10^{-6}$ years in oropouche S)
- **Cause of the narrow posteriors**: these trees hold chains of hundreds of zero-length branches, along which the backward messages lose accuracy and the posteriors shrink step by step ([H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md](H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md)). After the short-branch collapse of v0 ([M-timetree-pre-loop-short-branch-pruning-missing.md](M-timetree-pre-loop-short-branch-pruning-missing.md)), each build finishes at the same rate
- **A larger limit is no remedy**: at `--max-grid-points 2000000`, oropouche S stops at 2 037 322 points with half the spacing; at `30000000` it finishes after 838 seconds with a peak memory of 1.76 GB

## Rejected fixes

- **A point shortcut for a narrow branch, gated by the width of the node**: the shortcut picks the wrong mode of a node with two modes when one mode is narrower than the branch, and its gate has no error bound (measurement in the proposal, section "A point shortcut for a narrow branch picks the wrong mode of a two-mode node")
- **A product window cut where the product density falls below $10^{-15}$ of its peak**: a density cutoff gives no bound on the mass it removes. A distribution with neg-log ordinates $[0, 16 \ln 10]$ on $[0, 1]$ and a right tail slope of $10^{-20}$ holds almost all its mass in the tail; the cutoff at 0.9375 leaves 0.99999729 of the mass outside the window

## Fix direction

> [!IMPORTANT]
> **Decision required.** Moment matching of a narrow operand. An operand whose standard deviation is at most the grid spacing of the other operand is replaced by a shift to its mean and the three-point kernel $[s/2,\ 1 - s,\ s/2]$, $s = \sigma^2 / \Delta x^2$, where $\sigma$ is the standard deviation of the narrow operand and $\Delta x$ the spacing of the other operand. The kernel is exact for mass, mean and variance (proposal, section "Moment matching"), and the result keeps the spacing of the wide operand. This diverges from v0: v0 convolves a narrow branch by FFT at the spacing $\max(0.005\,m, \min(w_a, w_b)/150)$, which its floor $0.005\,m$ keeps coarse, and raises an error when the node support is more than 150 branch widths wide and the branch is narrower than 4 grid cells ([packages/legacy/treetime/treetime/node_interpolator.py#L162-L197](../../packages/legacy/treetime/treetime/node_interpolator.py#L162-L197); $m$ is the branch length of one mutation, $w_a$ and $w_b$ the widths at half maximum of the operands).

> [!IMPORTANT]
> **Decision required.** A product window from a bound on the lost mass. Before a product is evaluated at the finest spacing, bound the mass it would hold outside a candidate window, for example from the per-cell lower bound $\min y_f + \min y_g$ of the neg-log ordinates of the two factors over the cells of a coarse grid, whose exponentials bound the product mass of each cell, and grid only the window outside of which the bound is below $\epsilon$. Apply it only where the product is a forward posterior: windowing a backward message by its own mass removes mass that the posterior needs (proposal, section "Where to trim", and [kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md#a-hard-edge-at-an-fft-crop-point-removes-density-that-another-message-needs)).

## Validation

- The reproduction above finishes under the default grid point limit, every node dated
- `data/mpox/clade-ii/1000` with `--metadata-id-columns=accession --clock-rate=0.5 --max-iter=1` finishes, every node dated
- A forward convolution of a node posterior whose standard deviation is below the branch grid spacing with a branch distribution 100 years wide stays under the default grid point limit
- Dates and intervals against the current output on runs at realistic clock rates

## Related

- [kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md)
- [kb/decisions/distribution-grid-point-limit.md](../decisions/distribution-grid-point-limit.md)
- [M-timetree-near-flat-soft-tail-stretches-the-time-grid.md](M-timetree-near-flat-soft-tail-stretches-the-time-grid.md)
