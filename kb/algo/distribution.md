# Distribution and Convolution Algorithms

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

v1: `gaussian_convolution_pdf_grid()` (`#gaussian_convolution_pdf_grid`) in [`packages/treetime-analytical/src/gaussian.rs#L61-L63`](../../packages/treetime-analytical/src/gaussian.rs#L61-L63).

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

v1: `gaussian_product_params()` (`#gaussian_product_params`), `gaussian_product()` (`#gaussian_product`) in [`packages/treetime-analytical/src/gaussian.rs#L6-L50`](../../packages/treetime-analytical/src/gaussian.rs#L6-L50).

---

## Lazy Normalization Multiplication

Pointwise multiplication of discretized distributions with deferred normalization. It rescales when the running maximum drops below $10^{-100}$, reducing the number of normalization operations. Improved numerical accuracy has not been established. A non-finite or non-positive maximum currently produces an empty distribution, erasing the underlying numerical failure; the desired error policy is unresolved.

v1: `multiply_many_lazy_normalize()` (`#multiply_many_lazy_normalize`) in [`packages/treetime-ops/src/multiplication.rs#L39-L82`](../../packages/treetime-ops/src/multiplication.rs#L39-L82).

---

## Distribution Convolution

Polymorphic convolution dispatches across Point, Range, and Function variants. Function-Function convolution resamples both inputs to a common uniform grid and calls the FFT wrapper `treetime_ops::convolve_fft()`, which integrates each overlap by the trapezoid rule ([kb/decisions/distribution-convolution-trapezoid-rule.md](../decisions/distribution-convolution-trapezoid-rule.md)). The direct and ndarray-conv implementations in `treetime-ops` apply the same rule; distribution code does not use them. Point distributions short-circuit to a shift. Formula convolution returns an error.

v1: `distribution_convolution()` (`#distribution_convolution`) in [`packages/treetime-distribution/src/distribution_ops/convolve.rs#L13-L43`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L13-L43).

---

## Distribution Multiplication

Polymorphic pointwise multiplication dispatches across all `Distribution` variant pairs, including Formula. Range-Function, Function-Function, and Formula-Function multiplication use the exact analytical support intersection. Exact endpoint contact produces a Point distribution. Positive-width intersections use `round(overlap_width / dx) + 1` uniformly spaced points, at least 2; a count above the grid point limit stops the operation with an error ([kb/decisions/distribution-grid-point-limit.md](../decisions/distribution-grid-point-limit.md)); `dx` is the finer Function spacing for Function-Function and the concrete Function spacing for the other pairs. Formula-Formula retains a fixed 200-point grid because neither operand supplies a concrete spacing. v0 instead uses the union of original coordinates; [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md) records the intentional uniform-grid divergence. NaN and negative-roundoff semantics remain unresolved.

v1: `distribution_multiplication()` (`#distribution_multiplication`) in [`packages/treetime-distribution/src/distribution_ops/multiply.rs#L16`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L16).

## Edge Cell with an Offset Power Law

A grid edge with `BoundaryBehavior::HardEdgeCell` has zero density beyond the edge, and its outermost cell follows a power law whose origin lies at the offset $c$ beyond the edge $t_e$. With $d = |t - t_e| + c$ the distance from the origin, $d_a$ that of the inner cell end (the anchor) and $y_a$ the negative log density there, the cell holds $y(d) = y_a - b \ln(d/d_a) + s\,(d - d_a)$. The slope $s$ follows from the two stored cell ends, so the law reproduces both ends. The law stores the offset, not the origin, so a shift of the cell to absolute times keeps the distance exact near the edge. The input-mode Poisson branch likelihood uses it for its first cell ([kb/decisions/timetree-input-branch-likelihood-first-cell-and-bounds.md](../decisions/timetree-input-branch-likelihood-first-cell-and-bounds.md)).

The log density $b \ln d - s\,d$ is concave in $d$ for $b \ge 0$. It peaks at $d_m = b/s$ when $s > 0$ and this point lies inside the cell, which is then the mode of the cell; otherwise at a cell end. The mass of a part of the cell is integrated from the densest point $d_m$ of the part outward, in the variable $z = \ln(d/d_m)$, where the integrand $\exp\big((b+1)z - s\,d_m(e^z - 1)\big)$ is entire. Each panel is at most 0.5 wide, its log integrand changes by at most 1 across it, and its curvature times the squared width is at most 1; 12-point Gauss-Legendre quadrature then errs by about $10^{-14}$ of the panel mass. The integral stops where the density falls 40 nats below its peak: by concavity the part beyond holds less than $e^{-40}$ of the rest. The integrand is positive, so no sums cancel, and steep cells (slopes of $\pm 2 \cdot 10^6$ per unit distance) keep their relative precision. The inverse of the mass inside the cell uses Newton steps inside a bisection bracket; the inverse of the value bisects one monotone branch, between the mode and a cell end.

Products and quotients compose the laws of one side:

- **Same edge**: the exponents add (subtract for a quotient), and the offset is the exponent-weighted geometric mean, $\ln c^* = \sum_j b_j \ln c_j / \sum_j b_j$. With $a_j = \ln(c_j/c^*)$ and $\tau = c^*/(|t - t_e| + c^*)$, the true sum of log terms differs from the composed law by $\sum_j b_j \big[\ln(1 + (e^{a_j} - 1)\tau) - a_j \tau\big]$, because $\sum_j b_j a_j = 0$. Each bracket is the cumulant generating function of a Bernoulli variable minus its mean, at most $a_j^2/8$ by Hoeffding's lemma (<a id="cite-6"></a>Hoeffding 1963 [[6](#ref-6)]), so the error is at most $\sum_j b_j a_j^2 / 8$ nats, and 0 when the offsets agree. A quotient inverts a product exactly; a negative exponent is an internal error
- **Different edges**: the law at the edge of the result stays. The factor whose edge lies outside the result is linear between grid points in the result's first cell, with the interpolation error $b\,h^2 / (8 (\delta + c)^2)$ at a distance $\delta$ between the edges

v1: `struct EdgeCell` (`#EdgeCell`) in [`packages/treetime-grid/src/edge_cell.rs`](../../packages/treetime-grid/src/edge_cell.rs), `struct EdgeCellLaw` in [`packages/treetime-grid/src/edge_cell_law.rs`](../../packages/treetime-grid/src/edge_cell_law.rs), read by `struct MassProfile` in [`packages/treetime-distribution/src/distribution_ops/mass_profile.rs`](../../packages/treetime-distribution/src/distribution_ops/mass_profile.rs).

---

## Additional Algorithms

- Direct Convolution: `DirectConvolve`, `convolve_direct()` (`#convolve_direct`) in [`packages/treetime-ops/src/convolution.rs#L17-L39`](../../packages/treetime-ops/src/convolution.rs#L17-L39) - Explicit $O(nm)$ summation of the overlap sums, used by `treetime-validation` as an independent check of the FFT.
- Library Convolution: `NdarrayConvolve`, `convolve()` (`#convolve`) in [`packages/treetime-ops/src/convolution.rs#L41-L58`](../../packages/treetime-ops/src/convolution.rs#L41-L58) - Overlap sums from the ndarray-conv crate.
- Trapezoid rule: `trapezoid_from_overlap_sums()` (`#trapezoid_from_overlap_sums`) in [`packages/treetime-ops/src/convolution.rs#L99-L118`](../../packages/treetime-ops/src/convolution.rs#L99-L118) - Turns the overlap sums of every algorithm (direct, ndarray-conv, FFT) into trapezoid sums over each overlap; the two extreme outputs hold a quarter of the end products. Operands with fewer than 2 samples are rejected.
- Distribution Division: `distribution_division()` (`#distribution_division`) in [`packages/treetime-distribution/src/distribution_ops/divide.rs#L16`](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L16) - Cavity computation via pointwise division (subtraction in neg-log). Division reuses multiplication's support rule (`fn multiplication_support_intersection()`): per side the innermost hard bound when any operand is hard there, else the outermost soft bound. A soft divisor tail therefore extends the quotient -- sampled as bulk -- instead of truncating it, because the cavity dividend is a product containing the divisor as a factor so the quotient decays rather than spiking. Soft result sides are refit from the combined grid; a hard dividend edge bounds the quotient to zero (`Hard`), a divisor hard or `Error` edge yields `Error` beyond it. Exact endpoint contact is a Point distribution. [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md) defines the tails. Plain division floors small denominators at $10^{-10}$, changing the quotient; Formula division is unsupported (returns an error). The scientific small-divisor policy is unresolved.
- Distribution Negation: `distribution_negation()` (`#distribution_negation`) in [`packages/treetime-distribution/src/distribution_ops/negate.rs#L10-L18`](../../packages/treetime-distribution/src/distribution_ops/negate.rs#L10-L18) - Time reversal via `f(x) -> f(-x)`.

---

## Unimplemented

See [unimplemented](unimplemented.md) for full details:

- FFT Convolution with Delta Approximation
- Adaptive Simpson's Rule Convolution
- FWHM Computation

---

## References

- <a id="ref-1"></a>Cooley, James W., and John W. Tukey. 1965. "An Algorithm for the Machine Calculation of Complex Fourier Series." _Mathematics of Computation_ 19(90):297-301. https://doi.org/10.2307/2003354 [↩](#cite-1)
- <a id="ref-2"></a>Bromiley, Paul A. 2003. "Products and Convolutions of Gaussian Probability Density Functions." Tina Memo No. 2003-003. https://www.tina-vision.net/docs/memos/2003-003.pdf [↩](#cite-2)
- <a id="ref-3"></a>Ross, Sheldon M. 2014. _Introduction to Probability Models._ 11th ed. Academic Press. ISBN 978-0-12-407948-9. [↩](#cite-3)
- <a id="ref-4"></a>Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4(1):vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-4)
- <a id="ref-5"></a>Petersen, Kaare Brandt, and Michael Syskind Pedersen. 2012. _The Matrix Cookbook._ Technical report. https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf [↩](#cite-5)
- <a id="ref-6"></a>Hoeffding W. Probability inequalities for sums of bounded random variables. Journal of the American Statistical Association. 1963;58(301):13-30. doi:10.1080/01621459.1963.10500830 [↩](#cite-6)

---

## File Index

| File                                                                                                                   | Algorithms                                      |
| ---------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------- |
| [`packages/treetime-ops/src/convolution.rs`](../../packages/treetime-ops/src/convolution.rs)                           | Direct, ndarray-conv, FFT trapezoid convolution |
| [`packages/treetime-ops/src/multiplication.rs`](../../packages/treetime-ops/src/multiplication.rs)                     | Naive, aggressive, lazy multiplication          |
| [`packages/treetime-grid/src/edge_cell.rs`](../../packages/treetime-grid/src/edge_cell.rs) | Edge cell law: value, mode, mass by panel quadrature, inverses |
| [`packages/treetime-analytical/src/`](../../packages/treetime-analytical/src/)                                         | Gaussian, exponential analytical ops            |
| [`packages/treetime-distribution/src/distribution_core/`](../../packages/treetime-distribution/src/distribution_core/) | Distribution type                               |
| [`packages/treetime-distribution/src/distribution_ops/`](../../packages/treetime-distribution/src/distribution_ops/)   | Convolution, multiplication, division, negation |
