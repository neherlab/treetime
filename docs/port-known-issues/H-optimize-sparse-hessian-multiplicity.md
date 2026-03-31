# Sparse second derivative formula squares multiplicity incorrectly

The Newton-Raphson second derivative (Hessian) for sparse branch length optimization applies `.powi(2)` to the product `multiplicity * d1_term` instead of `d1_term` alone. For compressed site patterns with multiplicity $m$, the Hessian is inflated by factor $\sim m$, causing Newton to under-step on every edge.

## Location

[packages/treetime/src/commands/optimize/optimize_sparse_eval.rs#L36-L37](../../packages/treetime/src/commands/optimize/optimize_sparse_eval.rs#L36-L37)

```rust
second_derivative += multiplicity * (coefficients * &ev2_exp_ev).sum() / site_lh
  - (multiplicity * (coefficients * &ev_exp_ev).sum() / site_lh).powi(2);
```

## Mathematical background

Per-edge log-likelihood under the eigendecomposition formulation:

$$\ell(t) = \sum_i m_i \ln L_i(t), \quad L_i(t) = \sum_c k_{ic} e^{\lambda_c t}$$

where $m_i$ is the multiplicity of compressed site pattern $i$, $k_{ic}$ are precomputed coefficients from endpoint messages and GTR eigenvectors, and $\lambda_c$ are eigenvalues of the rate matrix $Q$.

First derivative:

$$\ell'(t) = \sum_i m_i \frac{\sum_c k_{ic} \lambda_c e^{\lambda_c t}}{L_i(t)}$$

Second derivative:

$$\ell''(t) = \sum_i m_i \left[ \frac{\sum_c k_{ic} \lambda_c^2 e^{\lambda_c t}}{L_i(t)} - \left(\frac{\sum_c k_{ic} \lambda_c e^{\lambda_c t}}{L_i(t)}\right)^2 \right]$$

The multiplicity $m_i$ is a linear factor on each site's contribution to $\ell''$. The squared term is the per-site first derivative squared, independent of multiplicity.

References:

- Felsenstein 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." J Mol Evol 17:368-376. https://doi.org/10.1007/BF01734359 -- original ML phylogenetic framework
- Stamatakis 2014. "RAxML Version 8." Bioinformatics 30(9):1312-1313. https://doi.org/10.1093/bioinformatics/btu033 -- eigendecomposition-based per-edge optimization with compressed site patterns

## Bug

The code computes:

$$m_i \cdot d_2 - (m_i \cdot d_1)^2 = m_i \cdot d_2 - m_i^2 \cdot d_1^2$$

instead of:

$$m_i \cdot (d_2 - d_1^2) = m_i \cdot d_2 - m_i \cdot d_1^2$$

where $d_1 = (\sum_c k_{ic} \lambda_c e^{\lambda_c t}) / L_i(t)$ and $d_2 = (\sum_c k_{ic} \lambda_c^2 e^{\lambda_c t}) / L_i(t)$.

The extra factor of $m_i$ on the squared term inflates $|\ell''|$. Since Newton step $\Delta t = -\ell'/\ell''$ divides by $\ell''$, a larger $|\ell''|$ produces smaller steps. Compressed fixed-site patterns in sparse representation carry multiplicities in the hundreds or thousands, so the inflation is large.

## Impact

All sparse branch length optimization is affected. Newton takes steps $\sim m$ times too small, requiring more iterations and converging to slightly wrong values if the iteration limit is reached. The dense path ([packages/treetime/src/commands/optimize/optimize_dense_eval.rs](../../packages/treetime/src/commands/optimize/optimize_dense_eval.rs)) does not have this bug because each site has implicit multiplicity 1.

This explains observed convergence differences between sparse and dense modes for the same dataset.

## Fix

Move `.powi(2)` inside so it applies only to the per-site derivative ratio:

```rust
second_derivative += multiplicity * (coefficients * &ev2_exp_ev).sum() / site_lh
  - multiplicity * ((coefficients * &ev_exp_ev).sum() / site_lh).powi(2);
```
