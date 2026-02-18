# Distribution and Convolution Algorithms

[Back to index](../index.md)

## FFT Convolution

| Property    | Value                                                                                                                  |
| ----------- | ---------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                             |
| v1 Location | `packages/treetime-ops/src/convolution.rs:33-37:`                                                                      |
| Reference   | Cooley & Tukey (1965). "An algorithm for the machine calculation of complex Fourier series." Math Comp, 19(90):297-301 |

O(n log n) via convolution theorem: `IFFT(FFT(f) * FFT(g))`.

---

## Gaussian Convolution

| Property    | Value                                                                                     |
| ----------- | ----------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                |
| v1 Location | `packages/treetime-analytical/src/gaussian.rs:104-117:`                                   |
| Reference   | Bromiley, P. (2003). "Products and Convolutions of Gaussian PDFs." Tina Memo No. 2003-003 |

Closed-form: `G1*G2 ~ N(mu1+mu2, sqrt(sigma1^2+sigma2^2))`.

---

## Exponential Convolution

| Property    | Value                                                            |
| ----------- | ---------------------------------------------------------------- |
| Type        | Well-known                                                       |
| v1 Location | `packages/treetime-analytical/src/exponential.rs:28-37:`         |
| Reference   | Ross, S. (2014). "Introduction to Probability Models." Chapter 5 |

Closed-form with special case for equal rates (Erlang-2).

---

## Gaussian-Exponential Convolution

| Property    | Value                                                             |
| ----------- | ----------------------------------------------------------------- |
| Type        | Well-known                                                        |
| v1 Location | `packages/treetime-analytical/src/gaussian_exponential.rs:12-14:` |
| Reference   | Sagulenko et al. (2018). TreeTime paper, Supplementary Section S2 |

Uses complementary error function (erfc).

---

## Gaussian Product

| Property    | Value                                                            |
| ----------- | ---------------------------------------------------------------- |
| Type        | Well-known                                                       |
| v1 Location | `packages/treetime-analytical/src/gaussian.rs:30-63:`            |
| Reference   | Petersen & Pedersen (2012). "The Matrix Cookbook." Section 8.1.8 |

Precision-weighted mean: `sigma* = 1/sqrt(sum(1/sigma_i^2))`, `mu* = sigma*^2 * sum(mu_i/sigma_i^2)`.

---

## ScaledDistribution

| Property    | Value                                                                        |
| ----------- | ---------------------------------------------------------------------------- |
| Type        | Well-known (log-sum-exp)                                                     |
| v1 Location | `packages/treetime-distribution/src/scaled_distribution.rs:13-91:`           |
| Reference   | Bishop, C.M. (2006). "Pattern Recognition and Machine Learning." Section 2.2 |

Decomposes `P(x) = exp(log_scale) * inner(x)` with `max(inner) = 1.0`.

---

## Lazy Normalization Multiplication

| Property    | Value                                                |
| ----------- | ---------------------------------------------------- |
| Type        | Custom                                               |
| v1 Location | `packages/treetime-ops/src/multiplication.rs:39-82:` |

Normalizes only when max drops below 1e-100, reducing rounding error vs aggressive normalization.

---

## Additional Algorithms

- **Riemann Sum Convolution** (`convolution.rs:9-19`): Direct O(n\*m)
- **Direct Library Convolution** (`convolution.rs:24-28`): ndarray-conv
- **Distribution Division** (`distribution_division.rs`): Cavity computation
- **Distribution Negation** (`distribution_negation.rs`): Time reversal
- **Quantile/Inverse CDF** (`distribution.rs:242-309`): Trapezoidal CDF + interpolation

---

## Unimplemented

See [unimplemented](../unimplemented/index.md) for full details:

- FFT Convolution with Delta Approximation
- Adaptive Simpson's Rule Convolution
- FWHM Computation
- Branch Length Interpolator (Input Mode)

---

## File Index

| File                                          | Algorithms                             |
| --------------------------------------------- | -------------------------------------- |
| `packages/treetime-ops/src/convolution.rs`    | Riemann, direct, FFT convolution       |
| `packages/treetime-ops/src/multiplication.rs` | Naive, aggressive, lazy multiplication |
| `packages/treetime-analytical/src/*.rs`       | Gaussian, exponential analytical ops   |
| `packages/treetime-distribution/src/*.rs`     | Distribution type system, operations   |
