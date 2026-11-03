# FFT convolution of node and branch distributions converges at first order at hard ends

## v0 location

`NodeInterpolator.convolve_fft()` (`#NodeInterpolator`, `#convolve_fft`) [packages/legacy/treetime/treetime/node_interpolator.py#L161](../../packages/legacy/treetime/treetime/node_interpolator.py#L161), the sum at [packages/legacy/treetime/treetime/node_interpolator.py#L214-L227](../../packages/legacy/treetime/treetime/node_interpolator.py#L214-L227)

## Erratum

v0 samples the node and the branch distribution on grids with spacing $dt$, multiplies their FFTs, and takes `-np.log(fft_res) - np.log(dt)`. This is the rectangle sum $dt \sum_i f_i\, g_{n-i}$ of the samples. Each overlap end gets the full weight $dt$ instead of the half weight of the trapezoid rule, so the result has an error of order $dt$ wherever the integrand is large at an overlap end. That is the case at a hard end with a large density: a branch without substitutions has its likelihood maximum at length 0, and a node distribution under a coalescent prior is pushed against the date of its earliest child. The error enters every message of the backward pass and the log marginal likelihood of node times, and it grows with the spacing, so at slow clock rates.

## Evidence

- **Contradicts the reference**: the convolution of the densities $f(t) = a e^{-at}$ and $g(t) = b e^{-bt}$ for $t \ge 0$ has the closed form $\frac{ab}{b-a}\left(e^{-at} - e^{-bt}\right)$. v0 `Distribution` objects of $f$ ($a = 1.5$, node) and $g$ ($b = 4$, branch) on $[0, 12]$ (4001 points, density 0 below 0, so hard left ends), convolved by `NodeInterpolator.convolve_fft()` with `inverse_time=True` and `one_mutation = 0`, so that $dt = \min(\mathrm{fwhm}) / \texttt{fft\_grid\_size}$. Error of the result density at $t$ and of its trapezoid mass against $(1 - e^{-12a})(1 - e^{-12b})$:

  | `fft_grid_size` | $dt$ | error at $t = 0.5$ | at $t = 1.0$ | at $t = 2.5$ | mass error |
  | --- | --- | --- | --- | --- | --- |
  | 25 | 6.960e-3 | 1.270e-2 | 5.052e-3 | 4.934e-4 | 1.914e-2 |
  | 50 | 3.480e-3 | 6.346e-3 | 2.523e-3 | 2.464e-4 | 9.570e-3 |
  | 100 | 1.740e-3 | 3.172e-3 | 1.261e-3 | 1.231e-4 | 4.785e-3 |
  | 200 | 8.700e-4 | 1.586e-3 | 6.303e-4 | 6.152e-5 | 2.393e-3 |
  | 400 | 4.350e-4 | 7.931e-4 | 3.151e-4 | 3.076e-5 | 1.196e-3 |
  | 800 | 2.175e-4 | 3.965e-4 | 1.576e-4 | 1.538e-5 | 5.981e-4 |

  The error halves with each halving of $dt$: the observed order between consecutive sizes is 1.000 to 1.002 in every column. A second-order quadrature reduces it four times per halving
- **Adjacent code is correct**: v0 integrates densities by the trapezoid rule in `Distribution.integrate_trapez()` (`#integrate_trapez`) [packages/legacy/treetime/treetime/distribution.py#L408](../../packages/legacy/treetime/treetime/distribution.py#L408), which gives each end sample the half weight that `convolve_fft()` leaves out

## v0 impact

The node-time messages and the log marginal likelihood of node times carry a first-order error at hard ends. Its size grows with the spacing, so with slow clock rates, and it biases the likelihood of node times toward slow clocks wherever the clock rate is fitted to it.

### Root date of a chain of zero-length branches

A caterpillar of $n$ tips with one sampling date, no substitution on any branch, and no coalescent prior has a closed-form root age: the second-largest of $n$ independent exponential ages of rate $\lambda$, with the CDF $(1 - u)^{n-1}(1 + (n-1)u)$, $u = e^{-\lambda x}$, and the mean $\sum_{m=2}^{n} 1/(\lambda m)$ ([kb/issues/H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md](../issues/H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md)). Every branch has its likelihood maximum at length 0, so every convolution of the backward pass has a large integrand at an overlap end. v0 0.12.1 (Nextstrain image of the pipelandia round, `treetime --sequence-length 754 --clock-rate 0.0014 --keep-root --keep-polytomies --time-marginal always`, $\lambda = 1.0556$ per year, tips at 2023.0) on a copy of `convolve_fft()` with $dt$ divided by 2 and 4 (both the FWHM term and the floor `one_mutation * 0.005`, which sets $dt$ in half of the convolutions), and with the end samples of each overlap sum at half weight. Error of the posterior mean of the root date in years, positive when the root is too young:

| $n$ | default $dt$ | $dt/2$ | $dt/4$ | half weight | half weight, $dt/2$ | half weight, $dt/4$ |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 50 | 0.1220 | 0.0618 | 0.0312 | 0.0042 | 0.00085 | 0.00020 |
| 200 | 1.2865 | 1.2494 | 1.2263 | 1.1501 | 1.1886 | 1.1796 |

At 50 tips the error halves with $dt$, and the half weight removes 97% of it and leaves an error that falls four to five times per halving: the rectangle sum is the cause. At 200 tips the rectangle sum causes 0.14 years of the 1.29 years, and the rest neither falls with $dt$ nor with the half weight, so a second v0 error dominates on long chains (see the issue above). The 90% max posterior region moves with the mean, and the standard deviation is within 0.006 years of the closed form in every run. Multiplying `fft_grid_size` by 4 without the floor gives 0.1123 at 50 tips, because the floor keeps $dt$ in half of the convolutions. v1 at its default grid has errors of 0.106 (50 tips) and 0.777 years (200 tips), and 0.0005 and 0.008 years at 16 times the points of the branch grid.

## v1 status

v1 integrates each overlap by the trapezoid rule (`fn trapezoid_from_overlap_sums()` in [packages/treetime-ops/src/convolution.rs](../../packages/treetime-ops/src/convolution.rs)) and converges at second order against the same closed form ([kb/decisions/distribution-convolution-trapezoid-rule.md](../decisions/distribution-convolution-trapezoid-rule.md)). Operands resampled to a common spacing keep their denser end on a grid point, which removes the remaining first-order error at the zero-length end of short branches ([kb/decisions/distribution-convolution-trapezoid-rule.md](../decisions/distribution-convolution-trapezoid-rule.md#operands-resampled-to-a-common-spacing-keep-their-denser-end-on-a-grid-point)).
