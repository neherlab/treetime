# Distribution and Convolution Algorithms

[Back to index](_index.md)

## FFT Convolution

| Property    | Value                                                                                                                                             |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                                        |
| v1 Location | [`packages/treetime-ops/src/convolution.rs#L33-L37`](../../packages/treetime-ops/src/convolution.rs#L33-L37) - `convolve_fft()` (`#convolve_fft`) |
| Reference   | Cooley & Tukey (1965). "An algorithm for the machine calculation of complex Fourier series." Math Comp, 19(90):297-301                            |

O(n log n) via convolution theorem: `IFFT(FFT(f) * FFT(g))`.

---

## Gaussian Convolution

| Property    | Value                                                                                                                                                                         |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                                                                    |
| v1 Location | [`packages/treetime-analytical/src/gaussian.rs#L110-L117`](../../packages/treetime-analytical/src/gaussian.rs#L110-L117) - `gaussian_convolution()` (`#gaussian_convolution`) |
| Reference   | Bromiley, P. (2003). "Products and Convolutions of Gaussian PDFs." Tina Memo No. 2003-003                                                                                     |

Closed-form: `G1*G2 ~ N(mu1+mu2, sqrt(sigma1^2+sigma2^2))`.

---

## Exponential Convolution

| Property    | Value                                                                                                                                                                                 |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                                                                            |
| v1 Location | [`packages/treetime-analytical/src/exponential.rs#L28-L37`](../../packages/treetime-analytical/src/exponential.rs#L28-L37) - `exponential_convolution()` (`#exponential_convolution`) |
| Reference   | Ross, S. (2014). "Introduction to Probability Models." Chapter 5                                                                                                                      |

Closed-form with special case for equal rates (Erlang-2).

---

## Gaussian-Exponential Convolution

| Property    | Value                                                                                                                                                                                                                     |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                                                                                                                |
| v1 Location | [`packages/treetime-analytical/src/gaussian_exponential.rs#L12-L14`](../../packages/treetime-analytical/src/gaussian_exponential.rs#L12-L14) - `gaussian_exponential_convolution()` (`#gaussian_exponential_convolution`) |
| Reference   | Sagulenko et al. (2018). TreeTime paper, Supplementary Section S2                                                                                                                                                         |

Uses complementary error function (erfc).

---

## Gaussian Product

| Property    | Value                                                                                                                                                                                                                       |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                                                                                                                  |
| v1 Location | [`packages/treetime-analytical/src/gaussian.rs#L30-L63`](../../packages/treetime-analytical/src/gaussian.rs#L30-L63) - `gaussian_product_params()` (`#gaussian_product_params`), `gaussian_product()` (`#gaussian_product`) |
| Reference   | Petersen & Pedersen (2012). "The Matrix Cookbook." Section 8.1.8                                                                                                                                                            |

Precision-weighted mean: `sigma* = 1/sqrt(sum(1/sigma_i^2))`, `mu* = sigma*^2 * sum(mu_i/sigma_i^2)`.

---

## ScaledDistribution

| Property    | Value                                                                                                                                                                                               |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known (log-sum-exp)                                                                                                                                                                            |
| v1 Location | [`packages/treetime-distribution/src/distribution_scaled/scaled.rs#L13`](../../packages/treetime-distribution/src/distribution_scaled/scaled.rs#L13) - `ScaledDistribution` (`#ScaledDistribution`) |
| Reference   | Bishop, C.M. (2006). "Pattern Recognition and Machine Learning." Section 2.2                                                                                                                        |

Decomposes `P(x) = exp(log_scale) * inner(x)` with `max(inner) = 1.0`.

---

## Lazy Normalization Multiplication

| Property    | Value                                                                                                                                                                                   |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Custom                                                                                                                                                                                  |
| v1 Location | [`packages/treetime-ops/src/multiplication.rs#L39-L82`](../../packages/treetime-ops/src/multiplication.rs#L39-L82) - `multiply_many_lazy_normalize()` (`#multiply_many_lazy_normalize`) |

Normalizes only when max drops below 1e-100, reducing rounding error vs aggressive normalization.

---

## Distribution Convolution

| Property    | Value                                                                                                                                                                                                                   |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Custom                                                                                                                                                                                                                  |
| v1 Location | [`packages/treetime-distribution/src/distribution_ops/convolve.rs#L13-L43`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L13-L43) - `distribution_convolution()` (`#distribution_convolution`) |

Polymorphic convolution dispatching across `Distribution` variants (Point, Range, Function). Resamples to uniform grid when inputs differ, delegates to `treetime_ops::convolve()` for Function-Function case.

---

## Distribution Multiplication

| Property    | Value                                                                                                                                                                                                                         |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Custom                                                                                                                                                                                                                        |
| v1 Location | [`packages/treetime-distribution/src/distribution_ops/multiply.rs#L14-L53`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L14-L53) - `distribution_multiplication()` (`#distribution_multiplication`) |

Polymorphic pointwise multiplication dispatching across all `Distribution` variant pairs including `Formula`. Formula distributions are discretized on the partner's grid before multiplication.

---

## Scaled Distribution Operations

| Property    | Value                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Custom                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| v1 Location | [`packages/treetime-distribution/src/distribution_scaled/multiply.rs#L51-L122`](../../packages/treetime-distribution/src/distribution_scaled/multiply.rs#L51-L122) - `scaled_distribution_multiplication()` (`#scaled_distribution_multiplication`), `scaled_distribution_multiply_many()` (`#scaled_distribution_multiply_many`); [`packages/treetime-distribution/src/distribution_scaled/convolve.rs#L9-L28`](../../packages/treetime-distribution/src/distribution_scaled/convolve.rs#L9-L28) - `scaled_distribution_convolution()` (`#scaled_distribution_convolution`) |

Log-scale-aware wrappers for multiplication and convolution on `ScaledDistribution`. Multi-way multiplication fast-paths through `treetime_ops::multiply_many_lazy_normalize()` when all inputs are grid-aligned Function distributions.

---

## Additional Algorithms

- **Riemann Sum Convolution** ([`packages/treetime-ops/src/convolution.rs#L9-L19`](../../packages/treetime-ops/src/convolution.rs#L9-L19)) - `convolve_riemann()` (`#convolve_riemann`): Direct O(n\*m)
- **Direct Library Convolution** ([`packages/treetime-ops/src/convolution.rs#L24-L28`](../../packages/treetime-ops/src/convolution.rs#L24-L28)) - `convolve()` (`#convolve`): ndarray-conv
- **Distribution Division** ([`packages/treetime-distribution/src/distribution_ops/divide.rs#L9-L51`](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L9-L51)) - `distribution_division()` (`#distribution_division`): Cavity computation via pointwise division
- **Scaled Distribution Division** ([`packages/treetime-distribution/src/distribution_scaled/divide.rs#L10-L33`](../../packages/treetime-distribution/src/distribution_scaled/divide.rs#L10-L33)) - `scaled_distribution_division()` (`#scaled_distribution_division`): Log-scale-aware wrapper
- **Distribution Negation** ([`packages/treetime-distribution/src/distribution_ops/negate.rs#L10-L18`](../../packages/treetime-distribution/src/distribution_ops/negate.rs#L10-L18)) - `distribution_negation()` (`#distribution_negation`): Time reversal via f(x) -> f(-x)
- **Quantile/Inverse CDF** ([`packages/treetime-distribution/src/distribution_core/distribution.rs#L236-L303`](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L236-L303)) - `quantile()` (`#quantile`): Trapezoidal CDF + linear interpolation
- **Confidence Interval** ([`packages/treetime-distribution/src/distribution_core/distribution.rs#L309-L313`](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L309-L313)) - `confidence_interval()` (`#confidence_interval`): Symmetric quantile bounds

---

## Unimplemented

See [unimplemented](unimplemented.md) for full details:

- FFT Convolution with Delta Approximation
- Adaptive Simpson's Rule Convolution
- FWHM Computation
- Branch Length Interpolator (Input Mode)

---

## File Index

| File                                                                                                                       | Algorithms                                          |
| -------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------- |
| [`packages/treetime-ops/src/convolution.rs`](../../packages/treetime-ops/src/convolution.rs)                               | Riemann, direct, FFT convolution                    |
| [`packages/treetime-ops/src/multiplication.rs`](../../packages/treetime-ops/src/multiplication.rs)                         | Naive, aggressive, lazy multiplication              |
| [`packages/treetime-analytical/src/`](../../packages/treetime-analytical/src/)                                             | Gaussian, exponential analytical ops                |
| [`packages/treetime-distribution/src/distribution_core/`](../../packages/treetime-distribution/src/distribution_core/)     | Distribution type, quantile, confidence interval    |
| [`packages/treetime-distribution/src/distribution_ops/`](../../packages/treetime-distribution/src/distribution_ops/)       | Convolution, multiplication, division, negation     |
| [`packages/treetime-distribution/src/distribution_scaled/`](../../packages/treetime-distribution/src/distribution_scaled/) | ScaledDistribution, scaled multiply/convolve/divide |
