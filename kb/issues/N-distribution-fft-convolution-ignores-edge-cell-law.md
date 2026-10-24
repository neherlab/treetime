# FFT convolution integrates an edge cell from its grid values only

> [!IMPORTANT]
> **Investigation required.** The change of timetree results from an exact edge cell in function-function convolution is not measured. Measure it on input-mode runs with short branches (`data/mpox/clade-ii/100`, `data/sc2/5100`): node dates against their 90% intervals and the log marginal likelihood of node times, with the exact cell mass in convolution against the grid-value quadrature.

A distribution with a `BoundaryBehavior::HardEdgeCell` side follows `EdgeCellLaw` in its outermost cell: a power law in the distance to an origin beyond the grid edge ([kb/decisions/timetree-input-branch-likelihood-first-cell-and-bounds.md](../decisions/timetree-input-branch-likelihood-first-cell-and-bounds.md)). The input-mode Poisson branch likelihood uses it for its first cell, where the cost $x/\mu - k \ln(x + c)$ curves on the scale $c = 10^{-3}\mu$, far below the grid step.

Function-function convolution ignores this law:

- `fn convolution_function_function_fine()` [packages/treetime-distribution/src/distribution_ops/convolve.rs#L264-L327](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L264-L327) resamples both operands to the finer spacing with `fn DistributionFunction::resample_dx()`, which interpolates through `fn GridFn::resample_range_dx_clamped()` [packages/treetime-grid/src/grid_fn.rs#L449-L460](../../packages/treetime-grid/src/grid_fn.rs#L449-L460) and `fn GridFn::interp()` without boundary laws, so the cell is linear between its grid values
- `fn convolve_fft()` then sums the grid samples, so the first cell of a branch contributes the quadrature of its two grid values instead of the mass of its law

In the time inference, every backward and forward convolution of a branch with a node distribution that is a function takes this path. The mass profile, the convolution window, and the convolutions with points and ranges read the law; the shift by the exact date of a leaf keeps it, so the products and posteriors at the parents of such leaves carry it.

## Size

For a branch with $k$ mutations, cell width $h$ (in units of $\mu$) and density $(x + c)^k e^{-x}$, the quadrature of the two grid values misses part of the first-cell mass when $0 < k < 1$. An estimate, not a measurement: at $k = 0.5$, $h = 0.1$ the quadrature gives 0.0160 against about 0.0199, 0.45% of the branch mass $\Gamma(1.5)$; at $k = 1$ the two agree within $10^{-4}$. At $k = 0$, which covers zero-length branches, the law is linear and nothing changes.

The measured effect of the channels that read the law is small: at most 4.9% of the 90% interval of a node date, and at most 1.94 days on the datasets measured ([kb/decisions/timetree-input-branch-likelihood-first-cell-and-bounds.md](../decisions/timetree-input-branch-likelihood-first-cell-and-bounds.md#effect-on-timetree-results)).

## Fix direction

Integrate the edge-cell mass of the law in the convolution quadrature once `dev` holds the trapezoid-rule rewrite of convolution quadrature from `rust` (`a1e98603f`, `dec0fe843`). A change to the current FFT path would conflict with that rewrite.
