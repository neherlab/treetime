# Clock lacks a rate uncertainty without covariation, and node marginals for the timetree start

Phase 0b of [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md).

## Problem

- `fn determine_rate_std()` ([`confidence.rs`](../../packages/treetime/src/timetree/confidence.rs)) gives a rate
  standard deviation only with `--covariation`; the augur node data report $\sqrt{(H^{-1})_{00}}$ of the unscaled
  regression Hessian as `rate_std` in every case
  ([`augur_node_data_refine.rs`](../../packages/app-output/src/augur_node_data_refine.rs)).
- The timetree loop needs node-time marginals to seed the expected lineage counts (phase 1).

## Fix direction

Gaussian message passing at the regression rate, following `gaussian_timetree.py`: a branch with $d$ mutations has
elapsed time $\mathcal N(\nu/\eta, \nu/\eta^2)$, $\nu = d + 1$, $\eta = \mu (D + M + N)/D$; messages $(S, R, W)$ with
mean $S - R/\eta$ and variance $W/\eta^2$; upward and outside passes give the node marginals, and the Hessian of
$\log L$ the rate and root covariance, converted to $\mu$. The regression's rate variance enters the node widths.

## Validation

Node RMSE, $\pm 1.645\sigma$ coverage and lineage-curve distance reproduce the proposal's evidence (240 validation
sets); `--confidence` interval coverage on the validation data.
