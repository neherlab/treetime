# Analytical Hessian loses precision via catastrophic cancellation

## Problem

The per-site second derivative of the branch-length log-likelihood computes variance via the difference-of-squares formula:

$$\ell''_i = \frac{\sum_c k_{ic} \lambda_c^2 e^{\lambda_c t}}{L_i} - \left(\frac{\sum_c k_{ic} \lambda_c e^{\lambda_c t}}{L_i}\right)^2$$

This is $E[\lambda^2] - E[\lambda]^2$ (variance of eigenvalues weighted by the site's coefficient distribution). When both terms are large and close, the subtraction loses significant digits. The result can have only a few correct digits despite each term being computed to full f64 precision.

Implementation: `packages/treetime/src/commands/optimize/optimize_eval.rs:41` in `evaluate_site_contributions()`:

```
second_derivative += m * (&coefficients * &ev2_exp_ev).sum() / site_lh - m * d1.powi(2);
```

where `d1 = (&coefficients * &ev_exp_ev).sum() / site_lh` is the first-derivative-to-likelihood ratio.

## Evidence

Property tests at `packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_prop_invariants.rs` require `max_relative = 1e-2` for the Hessian to agree with finite-difference approximation. The theoretical bound for central-difference second derivatives on smooth functions is $O(h^2) \sim 10^{-6}$ to $10^{-8}$. The 4-order-of-magnitude gap ($10^{-2}$ vs $10^{-6}$) indicates the analytical formula loses ~4 digits to cancellation.

## Impact

A 1% error in the Hessian produces a 1% error in the Newton step per iteration. Affects all three Newton methods (`newton`, `newton-sqrt`, `newton-log`). Does not affect Brent methods (derivative-free).

The $\ln(t)$ reparameterization does not fix this -- it fixes the indel Hessian conditioning but the substitution Hessian formula has this independent numerical issue.

## Stable alternative

Compute the numerator before dividing, avoiding the difference of two $O(1)$ quantities:

$$\ell''_i = \frac{\sum_c k_{ic} \lambda_c^2 e^{\lambda_c t} \cdot L_i - \left(\sum_c k_{ic} \lambda_c e^{\lambda_c t}\right)^2}{L_i^2}$$

The numerator $\sum k \lambda^2 e^{\lambda t} \cdot L_i - (\sum k \lambda e^{\lambda t})^2$ subtracts two products rather than two ratios, preserving more significant digits when $L_i$ is large.

An alternative approach: compute the variance using the one-pass Welford algorithm on the weighted eigenvalue distribution, but this changes the computation structure and may not vectorize as well with ndarray.

## Severity

Medium. Production bug affecting all Newton methods. The Hessian error is systematic (not random) and proportional to the condition number of the eigenvalue-coefficient distribution at each site.
