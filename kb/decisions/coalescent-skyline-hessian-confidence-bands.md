# Skyline confidence uses the full Hessian covariance

v1 reports curvature-based confidence bands for inferred coalescent time scales from the full penalized skyline Hessian. This differs from v0, which estimates each segment's curvature separately with finite differences.

**Type**: Scientific and implementation change from v0.

**Status**: Approved by the maintainer team.

**v0**: `MergerModel.optimize_skyline()` at [`packages/legacy/treetime/treetime/merger_models.py#L320-L333`](../../packages/legacy/treetime/treetime/merger_models.py#L320-L333) perturbs one $ln T_c$ value at a time by $0.1$. It uses the resulting diagonal finite-difference curvature and omits covariance caused by the skyline stiffness penalty.

**v1**: `optimize_skyline()` at [`packages/treetime/src/coalescent/skyline.rs`](../../packages/treetime/src/coalescent/skyline.rs) evaluates the exact symmetric tridiagonal Hessian at the optimum and solves against the identity matrix. The diagonal of the inverse gives the marginal variance for each $ln T_c$ value under a local <a id="gloss-use-1"></a>Laplace approximation <sup>[1](#gloss-1)</sup>. A normal approximation centered at one dominant mode uses the inverse observed curvature as covariance <a id="cite-1"></a>[Tierney and Kadane 1986](https://doi.org/10.1080/01621459.1986.10478240) [[1](#ref-1)]. Sparse neighboring dependence is represented by a <a id="gloss-use-2"></a>precision matrix <sup>[2](#gloss-2)</sup>, whose inverse is the covariance matrix <a id="cite-2"></a>[Rue et al. 2009](https://doi.org/10.1111/j.1467-9868.2008.00700.x) [[2](#ref-2)].

## Decision

At the optimized log time-scale vector, calculate the covariance and per-segment band as

$$
\Sigma = H^{-1}, \qquad \sigma_i = \sqrt{\Sigma_{ii}}, \qquad \left[T_{c,i} e^{-n_{\mathrm{std}}\sigma_i},\; T_{c,i} e^{n_{\mathrm{std}}\sigma_i}\right].
$$

where:

- $H$ -- Hessian of the penalized negative log-likelihood at the optimum
- $\Sigma$ -- local covariance matrix in $ln T_c$ coordinates
- $\sigma_i$ -- standard deviation of segment $i$ in $ln T_c$ coordinates
- $T_{c,i}$ -- inferred coalescent time scale for segment $i$
- $n_{\mathrm{std}}$ -- requested confidence-band width in standard deviations

Use the same calculation for a constant inferred $T_c$, represented as a one-segment skyline. A fixed user-supplied $T_c$ has no inferred variance and therefore no confidence band.

## Rationale

- The exact Hessian uses the derivatives already required by the Newton solve.
- The full inverse includes covariance from the adjacent-segment stiffness penalty. v0's one-coordinate finite differences cannot include this coupling.
- The multiplicative band is symmetric in $ln T_c$, preserves positive bounds, and matches v0's output form.

## Limits

The band is a local Gaussian approximation for the penalized model. It is not an exact confidence interval and does not establish repeated-sampling coverage. Its interpretation depends on the skyline likelihood, the Gaussian smoothing penalty, and the chosen $n_{\mathrm{std}}$.

## Glossary

1. <a id="gloss-1"></a> **Laplace approximation.** A local Gaussian approximation centered at a mode, with covariance from the inverse curvature at that mode. [↩](#gloss-use-1)
2. <a id="gloss-2"></a> **Precision matrix.** The inverse of a covariance matrix. Sparse precision entries encode conditional dependence between neighboring skyline segments. [↩](#gloss-use-2)

## References

1. <a id="ref-1"></a> Tierney, Luke, and Joseph B. Kadane. 1986. "Accurate Approximations for Posterior Moments and Marginal Densities." _Journal of the American Statistical Association_ 81:82-86. https://doi.org/10.1080/01621459.1986.10478240 [↩](#cite-1)
2. <a id="ref-2"></a> Rue, Håvard, Sara Martino, and Nicolas Chopin. 2009. "Approximate Bayesian Inference for Latent Gaussian Models by Using Integrated Nested Laplace Approximations." _Journal of the Royal Statistical Society Series B: Statistical Methodology_ 71:319-392. https://doi.org/10.1111/j.1467-9868.2008.00700.x [↩](#cite-2)
