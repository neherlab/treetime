# Distribution and Convolution Algorithms

[Back to index](README.md)

## FFT Convolution

O(n log n) convolution via the convolution theorem: `IFFT(FFT(f) * FFT(g))`. The discrete Fourier transform converts convolution (an O(n^2) operation in the time domain) into pointwise multiplication in the frequency domain, then transforms back.

v1: `convolve_fft()` (`#convolve_fft`) in [`packages/treetime-ops/src/convolution.rs#L33-L37`](../../packages/treetime-ops/src/convolution.rs#L33-L37).

Reference: Cooley & Tukey (1965). "An algorithm for the machine calculation of complex Fourier series." Math Comp, 19(90):297-301. doi:10.2307/2003354

---

## Gaussian Convolution

Closed-form convolution of two Gaussians: `G1 * G2 ~ N(mu1 + mu2, sqrt(sigma1^2 + sigma2^2))`. The convolution of two Gaussian PDFs is itself Gaussian with mean equal to the sum of means and variance equal to the sum of variances. No numerical integration required.

v1: `gaussian_convolution()` (`#gaussian_convolution`) in [`packages/treetime-analytical/src/gaussian.rs#L110-L117`](../../packages/treetime-analytical/src/gaussian.rs#L110-L117).

Reference: Bromiley (2003). "Products and Convolutions of Gaussian PDFs." Tina Memo No. 2003-003.

---

## Exponential Convolution

Closed-form convolution of two exponential distributions. When rates differ (a != b), the result is a hypoexponential distribution with PDF involving two exponential terms. When rates are equal (a == b), the result is the Erlang-2 distribution: `f(x) = a^2 * x * exp(-a*x)`.

v1: `exponential_convolution()` (`#exponential_convolution`) in [`packages/treetime-analytical/src/exponential.rs#L28-L37`](../../packages/treetime-analytical/src/exponential.rs#L28-L37).

Reference: Ross (2014). "Introduction to Probability Models." Chapter 5.

---

## Gaussian-Exponential Convolution

Convolution of a Gaussian with an exponential distribution, yielding a function involving the complementary error function (erfc). This arises in timetree inference when combining a Gaussian-like node time distribution with an exponentially distributed branch length.

v1: `gaussian_exponential_convolution()` (`#gaussian_exponential_convolution`) in [`packages/treetime-analytical/src/gaussian_exponential.rs#L12-L14`](../../packages/treetime-analytical/src/gaussian_exponential.rs#L12-L14).

Reference: Sagulenko, Puller & Neher (2018). "TreeTime." Virus Evolution, 4(1):vex042, Supplementary Section S2.

---

## Gaussian Product

Pointwise multiplication of two Gaussian PDFs (not convolution). The product of two Gaussians is an unnormalized Gaussian with precision-weighted mean:

```
sigma* = 1 / sqrt(1/sigma1^2 + 1/sigma2^2)
mu* = sigma*^2 * (mu1/sigma1^2 + mu2/sigma2^2)
```

This operation appears in belief propagation when combining independent messages at a node.

v1: `gaussian_product_params()` (`#gaussian_product_params`), `gaussian_product()` (`#gaussian_product`) in [`packages/treetime-analytical/src/gaussian.rs#L30-L63`](../../packages/treetime-analytical/src/gaussian.rs#L30-L63).

Reference: Petersen & Pedersen (2012). "The Matrix Cookbook." Section 8.1.8.

---

## ScaledDistribution

Decomposition `P(x) = exp(log_scale) * inner(x)` with `max(inner) = 1.0`. This log-sum-exp pattern prevents underflow when multiplying many small probabilities: the log-scale factor absorbs the magnitude while the inner distribution maintains full floating-point precision in the [0, 1] range. The standard technique for numerical stability in probabilistic computation (Bishop 2006).

v1: `ScaledDistribution` (`#ScaledDistribution`) in [`packages/treetime-distribution/src/distribution_scaled/scaled.rs#L13`](../../packages/treetime-distribution/src/distribution_scaled/scaled.rs#L13).

Reference: Bishop (2006). "Pattern Recognition and Machine Learning." Section 2.2.

---

## Lazy Normalization Multiplication

Pointwise multiplication of discretized distributions with deferred normalization. Normalizes only when the running maximum drops below 1e-100, reducing rounding error compared to aggressive per-step normalization. This is a TreeTime-specific optimization that trades a small number of underflow checks for fewer total division operations.

v1: `multiply_many_lazy_normalize()` (`#multiply_many_lazy_normalize`) in [`packages/treetime-ops/src/multiplication.rs#L39-L82`](../../packages/treetime-ops/src/multiplication.rs#L39-L82).

---

## Distribution Convolution

Polymorphic convolution dispatching across `Distribution` variants (Point, Range, Function, Formula). Resamples to a uniform grid when inputs have different grid spacings. Delegates to `treetime_ops::convolve()` for the Function-Function case (FFT or Riemann depending on size). Point distributions short-circuit to a simple shift operation.

v1: `distribution_convolution()` (`#distribution_convolution`) in [`packages/treetime-distribution/src/distribution_ops/convolve.rs#L13-L43`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L13-L43).

---

## Distribution Multiplication

Polymorphic pointwise multiplication dispatching across all `Distribution` variant pairs including `Formula`. Formula distributions (lazy closures) are discretized on the partner's grid before multiplication, producing a Function result. This is the mechanism by which coalescent contributions (represented as Formula) combine with discretized child messages during the backward pass.

v1: `distribution_multiplication()` (`#distribution_multiplication`) in [`packages/treetime-distribution/src/distribution_ops/multiply.rs#L14-L53`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L14-L53).

---

## Scaled Distribution Operations

Log-scale-aware wrappers for multiplication and convolution on `ScaledDistribution`. Multi-way multiplication fast-paths through `multiply_many_lazy_normalize()` when all inputs are grid-aligned Function distributions, avoiding repeated normalization overhead.

v1: `scaled_distribution_multiplication()` (`#scaled_distribution_multiplication`), `scaled_distribution_multiply_many()` (`#scaled_distribution_multiply_many`) in [`packages/treetime-distribution/src/distribution_scaled/multiply.rs#L51-L122`](../../packages/treetime-distribution/src/distribution_scaled/multiply.rs#L51-L122). `scaled_distribution_convolution()` (`#scaled_distribution_convolution`) in [`packages/treetime-distribution/src/distribution_scaled/convolve.rs#L9-L28`](../../packages/treetime-distribution/src/distribution_scaled/convolve.rs#L9-L28).

---

## Additional Algorithms

- **Riemann Sum Convolution**: `convolve_riemann()` (`#convolve_riemann`) in [`packages/treetime-ops/src/convolution.rs#L9-L19`](../../packages/treetime-ops/src/convolution.rs#L9-L19) - Direct O(n\*m) summation, used as reference implementation and for small arrays where FFT overhead dominates.
- **Direct Library Convolution**: `convolve()` (`#convolve`) in [`packages/treetime-ops/src/convolution.rs#L24-L28`](../../packages/treetime-ops/src/convolution.rs#L24-L28) - Delegates to ndarray-conv crate.
- **Distribution Division**: `distribution_division()` (`#distribution_division`) in [`packages/treetime-distribution/src/distribution_ops/divide.rs#L9-L51`](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L9-L51) - Cavity computation via pointwise division, used in the forward pass to compute the "outgroup message" (parent's information excluding the current child).
- **Scaled Distribution Division**: `scaled_distribution_division()` (`#scaled_distribution_division`) in [`packages/treetime-distribution/src/distribution_scaled/divide.rs#L10-L33`](../../packages/treetime-distribution/src/distribution_scaled/divide.rs#L10-L33) - Log-scale-aware wrapper.
- **Distribution Negation**: `distribution_negation()` (`#distribution_negation`) in [`packages/treetime-distribution/src/distribution_ops/negate.rs#L10-L18`](../../packages/treetime-distribution/src/distribution_ops/negate.rs#L10-L18) - Time reversal via `f(x) -> f(-x)`.
- **Quantile/Inverse CDF**: `quantile()` (`#quantile`) in [`packages/treetime-distribution/src/distribution_core/distribution.rs#L236-L303`](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L236-L303) - Trapezoidal CDF integration with linear interpolation for inverse lookup.
- **Confidence Interval**: `confidence_interval()` (`#confidence_interval`) in [`packages/treetime-distribution/src/distribution_core/distribution.rs#L309-L313`](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L309-L313) - Symmetric quantile bounds (e.g., 2.5% and 97.5% for a 95% CI).

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
