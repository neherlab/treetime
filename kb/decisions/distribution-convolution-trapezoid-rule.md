# Convolution of two time functions integrates each overlap by the trapezoid rule

**Type**: Intentional divergence from v0 (numerical method). v0 defect: [kb/v0-errata/timetree-fft-convolution-rectangle-rule.md](../v0-errata/timetree-fft-convolution-rectangle-rule.md).

The time passes convolve node-time distributions with branch-time distributions. v1 integrates each output of a convolution of two sampled functions by the trapezoid rule over the overlap of the operands. v0 takes the rectangle sum of the samples, which has a first-order error in the grid spacing wherever an operand has a large density at a grid end. That error biased the log marginal likelihood of node times (`log_lh_time`) at low signal and slow clock rates.

**v1 location**: `fn trapezoid_from_overlap_sums()` [packages/treetime-ops/src/convolution.rs#L99-L118](../../packages/treetime-ops/src/convolution.rs#L99-L118), applied to the overlap sums of every convolution algorithm of `treetime-ops`: `fn convolve_direct()` (`DirectConvolve`), `fn convolve()` (`NdarrayConvolve`) and `pub fn convolve_fft()` (`FftConvolve`) [packages/treetime-ops/src/convolution.rs#L17-L86](../../packages/treetime-ops/src/convolution.rs#L17-L86). The distribution code calls `convolve_fft()` from `fn convolution_function_function_fine()` [packages/treetime-distribution/src/distribution_ops/convolve.rs#L260](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L260).

**v0 location**: `NodeInterpolator.convolve_fft()` (`#convolve_fft`) [packages/legacy/treetime/treetime/node_interpolator.py#L161](../../packages/legacy/treetime/treetime/node_interpolator.py#L161), the sum at [packages/legacy/treetime/treetime/node_interpolator.py#L214-L227](../../packages/legacy/treetime/treetime/node_interpolator.py#L214-L227).

## What v1 computes

The operands are sampled on uniform grids with the same spacing:

$$
c_n \approx (f * g)(t_n) = \int f(s)\, g(t_n - s)\, ds, \qquad n = 0, \dots, F + G - 2
$$

where:

- $f_0, \dots, f_{F-1}$ and $g_0, \dots, g_{G-1}$: the samples of the two operands, $F, G \ge 2$
- $h$: the grid spacing
- $c_n$: the $n$-th output, at the grid point $t_n$
- $O_n = \{i_0, \dots, i_1\}$ with $i_0 = \max(0, n - G + 1)$ and $i_1 = \min(F - 1, n)$: the overlap of output $n$, the indices $i$ at which both $f_i$ and $g_{n-i}$ exist

Every output whose overlap has at least 2 samples ($1 \le n \le F + G - 3$) is the trapezoid rule over its overlap:

$$
c_n = h \sum_{i = i_0}^{i_1} f_i\, g_{n-i} - \frac{h}{2}\left(f_{i_0}\, g_{n-i_0} + f_{i_1}\, g_{n-i_1}\right)
$$

The first and the last output, whose overlap is one sample, hold $h f_0 g_0 / 4$ and $h f_{F-1} g_{G-1} / 4$ instead of the exact 0, because exact zeros fall below the trust floor of the tail reconstruction (`CONV_TRUST_FRACTION`, $10^{-13}$ of the peak, in `fn reconstruct_neg_log_tails()` [packages/treetime-distribution/src/distribution_ops/convolve.rs#L349](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L349)). Their effect on the total mass is second order in $h$. Operands with fewer than 2 samples are rejected with an error; a `DistributionFunction` always holds at least 2 samples, and a one-sample distribution is a `Point`, which the convolution shifts by exactly.

The rectangle sum of v0 is the first term alone:

$$
c_n^{\text{rect}} - c_n = \frac{h}{2}\left(f_{i_0}\, g_{n-i_0} + f_{i_1}\, g_{n-i_1}\right)
$$

By the Euler-Maclaurin formula, the trapezoid error of an integrand that is smooth up to the overlap ends is $O(h^2)$, so this difference is the $O(h)$ error of the rectangle sum ([Trefethen and Weideman 2014](https://doi.org/10.1137/130932132) [[1](#ref-1)]). It vanishes when the integrand is near 0 at both overlap ends and is large at a hard grid end with a large density.

## Why

In the backward pass both operands of most edge convolutions have a large density at a grid end:

- A branch without substitutions has its likelihood maximum at length 0
- The coalescent prior pushes a subtree distribution against the date of its earliest child

The grid spacing scales with $1/\mu$, where $\mu$ is the clock rate, so the rectangle error grows at slow rates. It biased `log_lh_time` (the log marginal likelihood $\log Z$ of node times) toward slow clocks. Under a coalescent prior the rate step moves the rate toward the maximum of $\log Z$, so the rate estimate was biased low ([kb/reports/timetree-low-signal-rate-bias-from-convolution-edges.md](../reports/timetree-low-signal-rate-bias-from-convolution-edges.md)).

## Evidence

- **Convergence order**: `test_convolution_function_function_hard_edges_converge_at_second_order` convolves $a e^{-at}$ and $b e^{-bt}$ ($a = 1.5$, $b = 4$, both on $[0, 4]$ with hard left ends) and compares the density with the closed form $\frac{ab}{b-a}\left(e^{-at} - e^{-bt}\right)$. `test_convolution_function_function_hard_edges_mass_converges_at_second_order` compares the trapezoid mass of the output with $(1 - e^{-4a})(1 - e^{-4b})$ (hard right ends too). Errors at 32, 64, 128 and 256 cells:

  | quantity | rectangle sum | trapezoid rule |
  | --- | --- | --- |
  | density at $t = 0.5$ | 2.345e-1, 1.156e-1, 5.738e-2, 2.859e-2 (order 1.00) | 6.572e-3, 1.645e-3, 4.114e-4, 1.029e-4 (order 2.00) |
  | density at $t = 1.0$ | 9.454e-2, 4.627e-2, 2.289e-2, 1.138e-2 (order 1.01) | 3.994e-3, 9.997e-4, 2.500e-4, 6.250e-5 (order 2.00) |
  | density at $t = 2.5$ | 9.294e-3, 4.533e-3, 2.238e-3, 1.112e-3 (order 1.01) | 4.577e-4, 1.146e-4, 2.865e-5, 7.163e-6 (order 2.00) |
  | total mass | 3.463e-1, 1.721e-1, 8.590e-2, 4.293e-2 (order 1.00) | 1.202e-2, 3.009e-3, 7.526e-4, 1.882e-4 (order 2.00) |

- **Hand values**: `test_convolution_algorithms_agree_with_the_trapezoid_rule_at_hard_ends` and `test_convolution_inner_outputs_are_trapezoid_sums_and_extreme_outputs_are_quarter_end_products` (`packages/treetime-ops/src/__tests__/test_convolution.rs`) check every algorithm against hand-computed trapezoid sums and the quarter end products
- **Analytic validation suite** (`treetime-validation`, every convolution algorithm, all within $3 \cdot 10^{-15}$ of each other). Suite `conv-exponential-exponential` (exponentials that start at the grid start): the mass error falls by 87 to 17000 times in 10 of 12 cases, for example from 0.015 to $1.2 \cdot 10^{-5}$ in `baseline_distinct_rates`, and the largest absolute error by 4 times in 11 cases. `tight_truncation` improves less (mass error 0.0083 to 0.0067), and `python_notebook_case` is unchanged, because its grid starts at $-1$ and the jump at 0 lies inside the grid. `slow_fast_decay` fails before and after for an unrelated reason ([kb/issues/N-validation-slow-fast-decay-expected-values-nan.md](../issues/N-validation-slow-fast-decay-expected-values-nan.md)). Suite `conv-gaussian-exponential`: 6 of 7 cases unchanged (the jump lies inside the grid). For the soft-end cases that become worse, see [Soft grid ends](#soft-grid-ends)
- **Exact posteriors of node times**: the reference of `test_marginal_exact_time_passes_match_the_integrated_posterior_without_a_prior` integrates the exact posterior by the trapezoid rule with the hard ends on grid points, converges at second order (observed order 2.000), and has an error of at most $1.6 \cdot 10^{-5}$. The largest CDF error of node R falls from $1.93 \cdot 10^{-3}$ to $7.2 \cdot 10^{-5}$ (case `mixed_counts`, limit $10^{-4}$) and from $1.42 \cdot 10^{-2}$ to $2.0 \cdot 10^{-4}$ (case `zero_counts_near_tips`, limit $10^{-3}$); both limits fail with the rectangle sum, and no case becomes worse. Remaining deviations: [kb/issues/M-timetree-time-pass-posteriors-deviate-from-the-exact-posterior.md](../issues/M-timetree-time-pass-posteriors-deviate-from-the-exact-posterior.md)
- **v0**: `NodeInterpolator.convolve_fft()` converges at order 1.000 to 1.002 against the same closed form ([kb/v0-errata/timetree-fft-convolution-rectangle-rule.md](../v0-errata/timetree-fft-convolution-rectangle-rule.md))
- **Log marginal likelihood of node times**: `test_marginal_likelihood_matches_the_reference_with_zero_count_branches_next_to_tips` (`packages/treetime/src/timetree/inference/__tests__/test_marginal_likelihood.rs`) computes $\log Z$ (the raw root log mass) of the tree `((A,B)I,C)root` with tips at 2015, 2014 and 2013 and Poisson branch likelihoods with $L = 1000$ and substitution counts A 1, B 0, I 1, C 0, so B and C are zero-count branches next to the tips. The oracle `reference_log_z()` (`packages/treetime/src/test_utils/reference_log_z.rs`) is an independent backward pass with per-cell quadrature and Richardson extrapolation over three spacings. It agrees with the closed forms of a Poisson cherry to about $10^{-14}$ and with `piecewise_log_z()` to $8 \cdot 10^{-9}$ or better, at an observed order of 1.99 to 2.00. Error $\lvert \log Z - \log Z_{\text{ref}} \rvert$:

  | prior | rate | trapezoid rule | rectangle sum |
  | --- | --- | --- | --- |
  | none | $10^{-4}$ | 3.0e-6 | 5.7e-4 |
  | none | $3 \cdot 10^{-4}$ | 2.7e-6 | 1.06e-3 |
  | none | $10^{-3}$ | 1.1e-6 | 1.18e-3 |
  | coalescent, $T_c = 5$ | $10^{-4}$ | 9.8e-6 | 2.98e-3 |
  | coalescent, $T_c = 5$ | $3 \cdot 10^{-4}$ | 4.2e-6 | 1.98e-3 |
  | coalescent, $T_c = 5$ | $10^{-3}$ | 1.4e-6 | 1.37e-3 |

  The trapezoid rule reduces the error 100 to 1000 times. The limits are $10^{-5}$ or $10^{-4}$ per case, and all 6 cases fail with the rectangle sum

## Operands resampled to a common spacing keep their denser end on a grid point

Both operands are resampled to the finer of their spacings before the FFT (`fn convolution_function_function_fine()` in [packages/treetime-distribution/src/distribution_ops/convolve.rs](../../packages/treetime-distribution/src/distribution_ops/convolve.rs)). When the range of an operand is not a whole number of cells of that spacing, one of its ends falls between grid points. `fn DistributionFunction::resample_dx()` keeps the end of higher density on a grid point, and the other end moves by less than a cell.

The grid used to be anchored at the left end. In the backward pass the negated branch distribution has its right end at branch length 0, where a short branch has its highest density. The convolution then lost (or, past half a cell, gained) the fraction $r$ of the edge cell, a relative mass error of $r\,h\,B(0) / \int B$ for the spacing $h$ and the branch density $B$: first order in $h$. On the tree of the $\log Z$ test above with B and I at zero substitutions, the I-to-root message missed the reference by $r = 0.020$ and $0.062$ of a cell at the rates $10^{-4}$ and $3 \cdot 10^{-4}$, and the predicted errors $3.37 \cdot 10^{-5}$ and $1.03 \cdot 10^{-4}$ match the measured $3.45 \cdot 10^{-5}$ and $1.04 \cdot 10^{-4}$. Every other step of that backward pass (leaf shifts, products, the root integral) was below $10^{-6}$.

Error of $\log Z$ against the reference over the 6 cases of the test, before and after anchoring at the denser end:

| substitution counts (A, B, I, C) | left anchor | denser-end anchor |
| --- | --- | --- |
| 1, 0, 1, 0 (one zero next to each of two tips) | 1.1e-6 to 9.8e-6 | unchanged |
| 0, 0, 1, 0 (all tip branches zero) | 3.8e-4 to 1.5e-3 | 1.5e-6 to 2.4e-5 |
| 1, 0, 0, 1 (B and I zero) | 3.2e-5 to 1.1e-4 | 3.1e-7 to 7.9e-6 |

With the denser-end anchor the I-to-root message converges at second order in the branch grid ($8.7 \cdot 10^{-7}$, $2.2 \cdot 10^{-7}$, $5.4 \cdot 10^{-8}$, $1.4 \cdot 10^{-8}$ at 6001 to 48001 points, rate $10^{-4}$). Tests: `test_convolution_keeps_the_dense_hard_end_of_a_resampled_operand_at_second_order` against the closed form of two exponentials, and the three count sets of `test_marginal_likelihood_matches_the_reference_with_zero_count_branches_next_to_tips`.

Other hard-edge effects measured against closed forms of two exponentials with hard ends at their maximum (aligned and misaligned grids, spacing ratios 1, 1.37 and 2, the backward layout through `fn convolve_across_edge()`):

- **Steep exponential edge**: an edge whose decay length is comparable to the spacing is resolved at second order once $kh < 1$ (log-decay rate $k = 40$: observed order 1.27 at $kh = 5$, 2.00 at $kh = 0.16$); it is the accepted second-order error of the grid ([discretization error of the default grid](distribution-tails-and-arithmetic.md#discretization-error-of-the-default-grid))
- **Extreme outputs** $h f_0 g_0 / 4$ at the point where two hard ends meet: the mass stays second order; the log density in the first cells next to that point is off by up to order 1, because a density that vanishes linearly is stored as log-linear pieces
- **`HardApproach` operands**: the convolution reads only grid values of the edge cell ([kb/issues/N-distribution-fft-convolution-ignores-edge-cell-law.md](../issues/N-distribution-fft-convolution-ignores-edge-cell-law.md))

## Soft grid ends

At a soft grid end (a `Linear` tail) the operand continues beyond the grid, and the convolution drops that tail. The extra half weight of the rectangle sum at the grid end partly offset the dropped mass, so Gaussians cut at a soft grid end are 1% to 26% worse in relative L2 error with the trapezoid rule: [kb/issues/M-distribution-convolution-drops-operand-tails.md](../issues/M-distribution-convolution-drops-operand-tails.md).

## Approval

The maintainers approved the divergence from v0 on the condition that the improvement is reproduced. The condition holds: the convergence tests, the validation suite, the exact-posterior tests and the $\log Z$ test above reproduce it, and the $\log Z$ test shows the error of the log marginal likelihood of node times fall 100 to 1000 times against an independent reference.

## Related

- [timetree-posterior-statistics-count-tail-mass.md](timetree-posterior-statistics-count-tail-mass.md): the mass model, whose trapezoid cells match the quadrature of the convolution
- [kb/issues/H-timetree-sparse-dense-node-times-differ.md](../issues/H-timetree-sparse-dense-node-times-differ.md): the effect of the trapezoid rule on the dense-to-sparse difference of node times
- [discretization error of the default grid](distribution-tails-and-arithmetic.md#discretization-error-of-the-default-grid): the log marginal likelihood at 300 grid points after this change

## References

1. <a id="ref-1"></a> Trefethen, Lloyd N., and J. A. C. Weideman. 2014. "The exponentially convergent trapezoidal rule." _SIAM Review_ 56 (3): 385-458. https://doi.org/10.1137/130932132
