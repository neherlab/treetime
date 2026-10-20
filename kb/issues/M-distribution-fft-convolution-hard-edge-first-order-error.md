# FFT convolution has a first-order error near a hard edge of an operand

The FFT convolution of a function with a function computes a Riemann sum of the sampled operands. When an operand jumps from zero to a finite density at a hard edge, the sum has a first-order error in the grid spacing $h$. v1 multiplies the discrete convolution of the sampled operands by the spacing in `fn convolve_fft()` [packages/treetime-ops/src/convolution.rs#L71-L84](../../packages/treetime-ops/src/convolution.rs#L71-L84), called by `fn convolution_function_function_fine()` [packages/treetime-distribution/src/distribution_ops/convolve.rs#L200](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L200). v0 has the same error, because it also takes a Riemann sum of the samples (`np.fft.irfft` of the sampled operands, then $-\ln dt$) [packages/legacy/treetime/treetime/node_interpolator.py#L214-L227](../../packages/legacy/treetime/treetime/node_interpolator.py#L214-L227).

> [!WARNING]
> **Needs review.** The measurement below comes from a stack of 2026-10-18 on top of `rust` of that day, with an exact mass model, an exact convolution of ranges, and the v0 grid rule of the FFT convolution, which `rust` does not have (last commit `38a325f2b`, local branch `worktree/fix/distribution-review-findings`). The test it describes exists only there. `rust` of 2026-10-19 still convolves by the Riemann sum above, so the first-order error at a hard edge applies to it; the numbers have not been measured on `rust`.

## Measurement

The test `test_convolution_exponential_gaussian_matches_closed_form_over_the_trusted_range` (`packages/treetime-distribution/src/distribution_ops/__tests__/test_convolve_function.rs`, lines 78 to 121, at commit `38a325f2b` only) convolves a unit Gaussian on ±8σ with the exponential $2e^{-2t}$ (hard left edge at 0, exact soft right tail), with the FFT spacing equal to the grid spacing $h$. It compares the neg-log result with the closed form $\tfrac12 a\,e^{-ax + a^2/2}\operatorname{erfc}((a - x)/\sqrt2)$ of `gaussian_exponential_convolution_grid()` at every stored sample whose closed-form density is above $10^{-13}$ of the peak. The closed form is the case of `test_scripts/conv_gauss_exp.ipynb`, where $a = 2$ is the rate of the exponential.

Far left of the bulk, the integrand is concentrated within about $1/(a + \lvert x\rvert)$ of the jump of the exponential, and the error grows to about $h(a + \lvert x\rvert)/2$ in neg-log. Largest error for $x < 10$, at the left trusted edge ($x \approx -7.2$ to $-7.4$):

| $h$   | error  |
| ----- | ------ |
| 0.04  | 0.182  |
| 0.02  | 0.0916 |
| 0.01  | 0.0458 |
| 0.005 | 0.0229 |

At the peak ($x \approx 0.43$) the error is 0.0050 at $h = 0.005$.

For $x \ge 10$ the integrand is far from the jump, and the result is exact up to FFT round-off: the density error stays below $10^{-15}$ of the peak density at every $h$. Relative to the small density near the crop level ($10^{-13}$ of the peak), this round-off is up to $2 \times 10^{-3}$ in neg-log.

On the stack, the test asserts first-order convergence of the left error and a bound of `1e-1` at the finest grid; the bound is looser than `1e-6` and cites this issue in a `TODO` line. A test of this kind on `rust` needs the closed form, `fn gaussian_exponential_convolution_grid()` of `treetime-analytical`, as a dev-dependency of `treetime-distribution`.

## Fix direction

> [!IMPORTANT]
> **Decision required.** Removing the error diverges from v0, which shares it. One option is to integrate the cell that contains the hard edge of an operand exactly instead of by the Riemann sum. The project default is v0 parity, so a fix needs user approval.
