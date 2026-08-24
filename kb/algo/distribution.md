# Distribution and Convolution Algorithms

[Back to index](README.md)

## FFT Convolution

O(n log n) convolution via the convolution theorem (<a id="cite-1"></a>[Cooley and Tukey 1965](https://doi.org/10.2307/2003354) [[1](#ref-1)]): `IFFT(FFT(f) * FFT(g))`. The discrete Fourier transform converts convolution (an O(n^2) operation in the time domain) into pointwise multiplication in the frequency domain, then transforms back.

v1: `convolve_fft()` (`#convolve_fft`) in [`packages/treetime-ops/src/convolution.rs#L33-L37`](../../packages/treetime-ops/src/convolution.rs#L33-L37).

---

## Gaussian Convolution

For $G_i(x)=A_i\exp[-(x-\mu_i)^2/(2\sigma_i^2)]$, closed-form Gaussian convolution gives (<a id="cite-2"></a>[Bromiley 2003](https://www.tina-vision.net/docs/memos/2003-003.pdf) [[2](#ref-2)])

$$\sigma_* = \sqrt{\sigma_1^2+\sigma_2^2},$$

$$
(G_1*G_2)(x)=A_1A_2\sqrt{2\pi}\frac{\sigma_1\sigma_2}{\sigma_*}
\exp\left[-\frac{(x-\mu_1-\mu_2)^2}{2\sigma_*^2}\right].
$$

For normalized Gaussian probability densities, this reduces to a normalized Gaussian whose mean and variance are the sums of the input means and variances.

v1: `gaussian_convolution()` (`#gaussian_convolution`) in [`packages/treetime-analytical/src/gaussian.rs#L124-L131`](../../packages/treetime-analytical/src/gaussian.rs#L124-L131).

---

## Exponential Convolution

Closed-form convolution of two exponential distributions (<a id="cite-3"></a>[Ross 2014](https://doi.org/10.1016/C2012-0-03564-8) [[3](#ref-3)]). When rates differ, the result is a hypoexponential distribution involving two exponential terms. When the rates differ by less than $10^{-15}$, the implementation evaluates the equal-rate Erlang-2 limit $f(x)=a^2x\exp(-ax)$.

v1: `exponential_convolution()` (`#exponential_convolution`) in [`packages/treetime-analytical/src/exponential.rs#L28-L37`](../../packages/treetime-analytical/src/exponential.rs#L28-L37).

---

## Gaussian-Exponential Convolution

The analytical crate implements convolution of a standard Gaussian and an exponential distribution with rate $a$:

$$h(x)=\frac{a}{2}\exp\left(-ax+\frac{a^2}{2}\right)\operatorname{erfc}\left(\frac{a-x}{\sqrt{2}}\right).$$

This operation is relevant to Gaussian and exponential message combinations in timetree theory (<a id="cite-4"></a>[Sagulenko, Puller, and Neher 2018](https://doi.org/10.1093/ve/vex042) [[4](#ref-4)], Supplementary Section S2), but current v1 timetree inference does not call this standalone analytical function directly.

v1: `gaussian_exponential_convolution()` (`#gaussian_exponential_convolution`) in [`packages/treetime-analytical/src/gaussian_exponential.rs#L12-L14`](../../packages/treetime-analytical/src/gaussian_exponential.rs#L12-L14).

---

## Gaussian Product

Pointwise multiplication of Gaussian factors produces an unnormalized Gaussian (<a id="cite-5"></a>[Petersen and Pedersen 2012](https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf) [[5](#ref-5)], Section 8.1.8). For factors with amplitudes $A_i$,

$$\sigma_*^{-2}=\sum_i\sigma_i^{-2}, \qquad \mu_*=\sigma_*^2\sum_i\frac{\mu_i}{\sigma_i^2},$$

$$\log A_* = \sum_i\log A_i - \frac12\sum_i\frac{(\mu_i-\mu_*)^2}{\sigma_i^2}.$$

This operation appears in belief propagation when combining independent messages at a node.

v1: `gaussian_product_params()` (`#gaussian_product_params`), `gaussian_product()` (`#gaussian_product`) in [`packages/treetime-analytical/src/gaussian.rs#L30-L63`](../../packages/treetime-analytical/src/gaussian.rs#L30-L63).

---

## ScaledDistribution

`ScaledDistribution` is a max-scaled representation $P(x)=e^s p(x)$. Normalizing constructors choose $s$ so that the sampled maximum of $p$ is one, which delays underflow when multiplying small values. `from_parts()` trusts the caller and does not enforce this invariant. This is max scaling, rather than a log-sum-exp calculation (<a id="cite-6"></a>[Bishop 2006](https://doi.org/10.1007/978-0-387-45528-0) [[6](#ref-6)], Section 2.2).

v1: `ScaledDistribution` (`#ScaledDistribution`) in [`packages/treetime-distribution/src/distribution_scaled/scaled.rs#L13`](../../packages/treetime-distribution/src/distribution_scaled/scaled.rs#L13).

---

## Lazy Normalization Multiplication

Pointwise multiplication of discretized distributions with deferred normalization. It rescales when the running maximum drops below $10^{-100}$, reducing the number of normalization operations. Improved numerical accuracy has not been established. A non-finite or non-positive maximum currently produces an empty distribution, erasing the underlying numerical failure; the desired error policy is unresolved.

v1: `multiply_many_lazy_normalize()` (`#multiply_many_lazy_normalize`) in [`packages/treetime-ops/src/multiplication.rs#L39-L82`](../../packages/treetime-ops/src/multiplication.rs#L39-L82).

---

## Distribution Convolution

Polymorphic convolution dispatches across Point, Range, and Function variants. Function-Function convolution resamples both inputs to a common uniform grid and calls the direct `ndarray-conv` wrapper `treetime_ops::convolve()`. The standalone FFT and Riemann implementations are not selected automatically. Point distributions short-circuit to a shift. Formula convolution currently panics rather than returning a typed error.

v1: `distribution_convolution()` (`#distribution_convolution`) in [`packages/treetime-distribution/src/distribution_ops/convolve.rs#L13-L43`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L13-L43).

---

## Distribution Multiplication

Polymorphic pointwise multiplication dispatches across all `Distribution` variant pairs, including Formula. Range-Function, Function-Function, and Formula-Function multiplication use the exact analytical support intersection. Exact endpoint contact produces a Point distribution. Positive-width intersections use `round(overlap_width / dx) + 1` uniformly spaced points, clamped to `[2, 1_000_000]`; `dx` is the finer Function spacing for Function-Function and the concrete Function spacing for the other pairs. Formula-Formula retains a fixed 200-point grid because neither operand supplies a concrete spacing. v0 instead uses the union of original coordinates; [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md) records the intentional uniform-grid divergence. NaN and negative-roundoff semantics remain unresolved.

v1: `distribution_multiplication()` (`#distribution_multiplication`) in [`packages/treetime-distribution/src/distribution_ops/multiply.rs#L16`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L16).

---

## Scaled Distribution Operations

Log-scale-aware wrappers for multiplication and convolution on `ScaledDistribution`. Multi-way multiplication fast-paths through `multiply_many_lazy_normalize()` when all inputs are grid-aligned Function distributions, avoiding repeated normalization overhead.

v1: `scaled_distribution_multiplication()` (`#scaled_distribution_multiplication`), `scaled_distribution_multiply_many()` (`#scaled_distribution_multiply_many`) in [`packages/treetime-distribution/src/distribution_scaled/multiply.rs#L51-L122`](../../packages/treetime-distribution/src/distribution_scaled/multiply.rs#L51-L122). `scaled_distribution_convolution()` (`#scaled_distribution_convolution`) in [`packages/treetime-distribution/src/distribution_scaled/convolve.rs#L9-L28`](../../packages/treetime-distribution/src/distribution_scaled/convolve.rs#L9-L28).

---

## Additional Algorithms

- Riemann Sum Convolution: `convolve_riemann()` (`#convolve_riemann`) in [`packages/treetime-ops/src/convolution.rs#L9-L19`](../../packages/treetime-ops/src/convolution.rs#L9-L19) - Explicit $O(nm)$ summation used as an independent validation oracle and available operation.
- Direct Library Convolution: `convolve()` (`#convolve`) in [`packages/treetime-ops/src/convolution.rs#L24-L28`](../../packages/treetime-ops/src/convolution.rs#L24-L28) - Delegates to ndarray-conv crate.
- Distribution Division: `distribution_division()` (`#distribution_division`) in [`packages/treetime-distribution/src/distribution_ops/divide.rs#L16`](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L16) - Cavity computation via pointwise division (subtraction in neg-log). Division reuses multiplication's support rule (`fn multiplication_support_intersection()`): per side the innermost hard bound when any operand is hard there, else the outermost soft bound. A soft divisor tail therefore extends the quotient -- sampled as bulk -- instead of truncating it, because the cavity dividend is a product containing the divisor as a factor so the quotient decays rather than spiking. Soft result sides are refit from the combined grid; a hard dividend edge bounds the quotient to zero (`Hard`), a divisor hard or `Error` edge yields `Error` beyond it. Exact endpoint contact is a Point distribution. [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md) defines the tails. Plain division floors small denominators at $10^{-10}$, changing the quotient; Formula division is unsupported (returns an error). The scientific small-divisor policy is unresolved.
- Scaled Distribution Division: `scaled_distribution_division()` (`#scaled_distribution_division`) in [`packages/treetime-distribution/src/distribution_scaled/divide.rs#L10-L33`](../../packages/treetime-distribution/src/distribution_scaled/divide.rs#L10-L33) - Log-scale-aware wrapper.
- Distribution Negation: `distribution_negation()` (`#distribution_negation`) in [`packages/treetime-distribution/src/distribution_ops/negate.rs#L10-L18`](../../packages/treetime-distribution/src/distribution_ops/negate.rs#L10-L18) - Time reversal via `f(x) -> f(-x)`.
- Quantile/Inverse CDF: `quantile()` (`#quantile`) in [`packages/treetime-distribution/src/distribution_core/distribution.rs#L236-L293`](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L236-L293) - Trapezoidal CDF integration with linear interpolation. Formula inputs and zero or non-finite integrated mass fall back to `likely_time()` instead of computing a quantile.
- Confidence Interval: `confidence_interval()` (`#confidence_interval`) in [`packages/treetime-distribution/src/distribution_core/distribution.rs#L295-L302`](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L295-L302) - Returns quantiles at the two caller-supplied probabilities; the interval is symmetric only when those probabilities are symmetric.

---

## Unimplemented

See [unimplemented](unimplemented.md) for full details:

- FFT Convolution with Delta Approximation
- Adaptive Simpson's Rule Convolution
- FWHM Computation
- Branch Length Interpolator (Input Mode)

---

## References

- <a id="ref-1"></a>Cooley, James W., and John W. Tukey. 1965. "An Algorithm for the Machine Calculation of Complex Fourier Series." _Mathematics of Computation_ 19(90):297-301. https://doi.org/10.2307/2003354 [↩](#cite-1)
- <a id="ref-2"></a>Bromiley, Paul A. 2003. "Products and Convolutions of Gaussian Probability Density Functions." Tina Memo No. 2003-003. https://www.tina-vision.net/docs/memos/2003-003.pdf [↩](#cite-2)
- <a id="ref-3"></a>Ross, Sheldon M. 2014. _Introduction to Probability Models._ 11th ed. Academic Press. ISBN 978-0-12-407948-9. [↩](#cite-3)
- <a id="ref-4"></a>Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4(1):vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-4)
- <a id="ref-5"></a>Petersen, Kaare Brandt, and Michael Syskind Pedersen. 2012. _The Matrix Cookbook._ Technical report. https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf [↩](#cite-5)
- <a id="ref-6"></a>Bishop, Christopher M. 2006. _Pattern Recognition and Machine Learning._ Springer. ISBN 978-0-387-31073-2. [↩](#cite-6)

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
