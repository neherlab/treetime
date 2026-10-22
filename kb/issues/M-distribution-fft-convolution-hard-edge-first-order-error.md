# FFT convolution has a first-order error near a hard edge of an operand

The convolution of two time functions integrates each overlap by the trapezoid rule ([kb/decisions/distribution-convolution-trapezoid-rule.md](../decisions/distribution-convolution-trapezoid-rule.md)), which removes the first-order error of the rectangle sum for operands with a large, smooth density at a grid end. The trapezoid rule is second order only for an integrand that is smooth on the grid scale up to the overlap end. Three cases stay worse. Code: `fn trapezoid_from_overlap_sums()` in [packages/treetime-ops/src/convolution.rs](../../packages/treetime-ops/src/convolution.rs), called through `convolve_fft()` by `fn convolution_function_function_fine()` in [packages/treetime-distribution/src/distribution_ops/convolve.rs](../../packages/treetime-distribution/src/distribution_ops/convolve.rs).

## Remaining cases

- **`HardApproach` operand**: a boundary whose density approaches the hard end by a power law with exponent $b$ inside the edge cell (`GridFn::hard_cell()` in [packages/treetime-grid/src/grid_fn.rs](../../packages/treetime-grid/src/grid_fn.rs)). The convolution reads only the grid samples, and the trapezoid rule over a cell with that law has an error of order $h^{1+b}$, so for $b < 1$ the convergence is below second order, where $h$ is the grid spacing
- **Steep exponential edge**: an operand whose log density falls with slope $k$ at a grid end. The per-cell trapezoid error is about $(kh)^2/12$, which is no longer small when $kh \gtrsim 1$. The slope of the coalescent prior factor has a fixed time scale while $h$ scales with $1/\mu$ (the clock rate $\mu$), so this error grows again at very slow rates
- **Extreme outputs**: the first and the last output, whose overlap is one sample, hold $h f_0 g_0/4$ and $h f_{F-1} g_{G-1}/4$ instead of the exact 0, where $f_i$ and $g_j$ are the samples of the two operands with $F$ and $G$ samples. Their effect on the total mass is second order in $h$, but for two operands with hard left ends the output mass changes by about $-h^2 f_0 g_0 / 8$, the same order as the trapezoid error. Exact zeros cannot be stored by the current code: `fn reconstruct_neg_log_tails()` crops outputs below `CONV_TRUST_FRACTION` ($10^{-13}$ of the peak), so a zero at a hard end would be cropped, the hard end of the result would move one cell inward, and a short result (4 outputs) would keep 2 trusted points and fail with "Convolution left too few trusted points to reconstruct tails", because its tail margin `CONV_TAIL_MARGIN.min(n / 3)` is 0

## Measurement where hard ends meet

`test_marginal_likelihood_matches_the_reference_with_zero_count_branches_next_to_tips` (`packages/treetime/src/timetree/inference/__tests__/test_marginal_likelihood.rs`) compares $\log Z$ of the tree `((A,B)I,C)root` (tips at 2015, 2014, 2013, Poisson branch likelihoods with $L = 1000$) with the independent oracle `reference_log_z()`, at the rates $10^{-4}$, $3 \cdot 10^{-4}$ and $10^{-3}$, without a prior and with a coalescent prior of $T_c = 5$. With one zero-count branch next to each of two tips (counts A 1, B 0, I 1, C 0) the error is at most $9.8 \cdot 10^{-6}$. Where several hard ends meet, it stays much larger, $\lvert \log Z - \log Z_{\text{ref}} \rvert$ over the same 6 cases:

| substitution counts (A, B, I, C) | trapezoid rule | rectangle sum |
| --- | --- | --- |
| 0, 0, 1, 0 (all tip branches zero) | 3.8e-4 to 1.5e-3 | 1.3e-3 to 4.9e-3 |
| 1, 0, 0, 1 (B and I zero) | 3.2e-5 to 1.1e-4 | 7.4e-4 to 2.5e-3 |

The `TODO` marker of the test points to this issue.

## Fix direction

- Integrate the edge cell exactly with its cell law (`HardApproach` power law, or the exponential of a steep edge) in an $O(N)$ end correction after the FFT, where $N$ is the number of outputs
- Or apply Gregory end corrections: fixed weights on a few samples near each overlap end that cancel further terms of the Euler-Maclaurin expansion of the trapezoid error ([Fornberg 2021](https://doi.org/10.1137/18M1229353))
- For the extreme outputs: store a zero at a hard end as the hard end of the result (an infinite negative log at the grid end) instead of a sample that the tail crop reads, and handle short results without the tail reconstruction

## Validation

Measure the convergence order of `distribution_convolution()` with a `HardApproach` operand ($b \in \{0.5, 1, 2\}$) and with an exponential edge at $kh \in \{0.1, 1, 5\}$ against closed forms, in the manner of `test_convolution_function_function_hard_edges_converge_at_second_order`. The fix holds when every case converges at second order, the extreme outputs hold the exact trapezoid value 0 without a failure on short results, and the two cases of the table above reach the error of the case with one zero-count branch per tip.
