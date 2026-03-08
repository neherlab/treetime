# Distribution and Convolution Algorithms

[Back to index](../index.md)

## FFT Convolution

| Property    | Value                                                                                                                                                |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                                           |
| v1 Location | [`packages/treetime-ops/src/convolution.rs#L33-L37`](../../../packages/treetime-ops/src/convolution.rs#L33-L37) - `convolve_fft()` (`#convolve_fft`) |
| Reference   | Cooley & Tukey (1965). "An algorithm for the machine calculation of complex Fourier series." Math Comp, 19(90):297-301                               |

O(n log n) via convolution theorem: `IFFT(FFT(f) * FFT(g))`.

---

## Gaussian Convolution

| Property    | Value                                                                                                                                                                            |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                                                                       |
| v1 Location | [`packages/treetime-analytical/src/gaussian.rs#L104-L117`](../../../packages/treetime-analytical/src/gaussian.rs#L104-L117) - `gaussian_convolution()` (`#gaussian_convolution`) |
| Reference   | Bromiley, P. (2003). "Products and Convolutions of Gaussian PDFs." Tina Memo No. 2003-003                                                                                        |

Closed-form: `G1*G2 ~ N(mu1+mu2, sqrt(sigma1^2+sigma2^2))`.

---

## Exponential Convolution

| Property    | Value                                                                                                                                                                                    |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                                                                               |
| v1 Location | [`packages/treetime-analytical/src/exponential.rs#L28-L37`](../../../packages/treetime-analytical/src/exponential.rs#L28-L37) - `exponential_convolution()` (`#exponential_convolution`) |
| Reference   | Ross, S. (2014). "Introduction to Probability Models." Chapter 5                                                                                                                         |

Closed-form with special case for equal rates (Erlang-2).

---

## Gaussian-Exponential Convolution

| Property    | Value                                                                                                                                                                                                                        |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                                                                                                                   |
| v1 Location | [`packages/treetime-analytical/src/gaussian_exponential.rs#L12-L14`](../../../packages/treetime-analytical/src/gaussian_exponential.rs#L12-L14) - `gaussian_exponential_convolution()` (`#gaussian_exponential_convolution`) |
| Reference   | Sagulenko et al. (2018). TreeTime paper, Supplementary Section S2                                                                                                                                                            |

Uses complementary error function (erfc).

---

## Gaussian Product

| Property    | Value                                                                                                                                                                                                                          |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Type        | Well-known                                                                                                                                                                                                                     |
| v1 Location | [`packages/treetime-analytical/src/gaussian.rs#L30-L63`](../../../packages/treetime-analytical/src/gaussian.rs#L30-L63) - `gaussian_product_params()` (`#gaussian_product_params`), `gaussian_product()` (`#gaussian_product`) |
| Reference   | Petersen & Pedersen (2012). "The Matrix Cookbook." Section 8.1.8                                                                                                                                                               |

Precision-weighted mean: `sigma* = 1/sqrt(sum(1/sigma_i^2))`, `mu* = sigma*^2 * sum(mu_i/sigma_i^2)`.

---

## ScaledDistribution

| Property    | Value                                                                                                                                                                                                                            |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known (log-sum-exp)                                                                                                                                                                                                         |
| v1 Location | [`packages/treetime-distribution/src/distribution_scaled/distribution_scaled.rs#L13`](../../../packages/treetime-distribution/src/distribution_scaled/distribution_scaled.rs#L13) - `ScaledDistribution` (`#ScaledDistribution`) |
| Reference   | Bishop, C.M. (2006). "Pattern Recognition and Machine Learning." Section 2.2                                                                                                                                                     |

Decomposes `P(x) = exp(log_scale) * inner(x)` with `max(inner) = 1.0`.

---

## Lazy Normalization Multiplication

| Property    | Value                                                                                                                                                                                      |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Type        | Custom                                                                                                                                                                                     |
| v1 Location | [`packages/treetime-ops/src/multiplication.rs#L39-L82`](../../../packages/treetime-ops/src/multiplication.rs#L39-L82) - `multiply_many_lazy_normalize()` (`#multiply_many_lazy_normalize`) |

Normalizes only when max drops below 1e-100, reducing rounding error vs aggressive normalization.

---

## Additional Algorithms

- **Riemann Sum Convolution** ([`packages/treetime-ops/src/convolution.rs#L9-L19`](../../../packages/treetime-ops/src/convolution.rs#L9-L19)) - `convolve_riemann()` (`#convolve_riemann`): Direct O(n\*m)
- **Direct Library Convolution** ([`packages/treetime-ops/src/convolution.rs#L24-L28`](../../../packages/treetime-ops/src/convolution.rs#L24-L28)) - `convolve()` (`#convolve`): ndarray-conv
- **Distribution Division** ([`packages/treetime-distribution/src/distribution_ops/divide.rs`](../../../packages/treetime-distribution/src/distribution_ops/divide.rs)): Cavity computation
- **Distribution Negation** ([`packages/treetime-distribution/src/distribution_ops/negate.rs`](../../../packages/treetime-distribution/src/distribution_ops/negate.rs)): Time reversal
- **Quantile/Inverse CDF** ([`packages/treetime-distribution/src/distribution_core/distribution.rs#L242-L309`](../../../packages/treetime-distribution/src/distribution_core/distribution.rs#L242-L309)) - `quantile()` (`#quantile`): Trapezoidal CDF + interpolation

---

## Unimplemented

See [unimplemented](../unimplemented/index.md) for full details:

- FFT Convolution with Delta Approximation
- Adaptive Simpson's Rule Convolution
- FWHM Computation
- Branch Length Interpolator (Input Mode)

---

## File Index

| File                                                                                                  | Algorithms                             |
| ----------------------------------------------------------------------------------------------------- | -------------------------------------- |
| [`packages/treetime-ops/src/convolution.rs`](../../../packages/treetime-ops/src/convolution.rs)       | Riemann, direct, FFT convolution       |
| [`packages/treetime-ops/src/multiplication.rs`](../../../packages/treetime-ops/src/multiplication.rs) | Naive, aggressive, lazy multiplication |
| [`packages/treetime-analytical/src/`](../../../packages/treetime-analytical/src/)                     | Gaussian, exponential analytical ops   |
| [`packages/treetime-distribution/src/`](../../../packages/treetime-distribution/src/)                 | Distribution type system, operations   |
