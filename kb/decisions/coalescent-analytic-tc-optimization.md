# Optimal constant Tc computed analytically instead of numerically

v1 computes the maximum-likelihood constant coalescent time scale $T_c$ from a
closed-form expression. v0 finds it by numerically minimizing the negative total
log-likelihood with Brent's method.

**Type**: Implementation change (same estimand, exact instead of iterative).

**v0 location**: `MergerModel.optimize_Tc()` at
[`packages/legacy/treetime/treetime/merger_models.py#L259-L279`](../../packages/legacy/treetime/treetime/merger_models.py). Runs
`scipy.optimize.minimize_scalar(cost, bracket=[-20.0, 2.0], method='brent')` where
`cost(logTc) = -total_LH()` re-evaluates the whole-tree likelihood for each
candidate $T_c$.

**v1 location**: `optimize_tc()` at
[`packages/treetime/src/coalescent/optimize_tc.rs`](../../packages/treetime/src/coalescent/optimize_tc.rs). The
$T_c$-independent integral $I$ reuses `compute_integral_merger_rate()` at
[`packages/treetime/src/coalescent/integration.rs`](../../packages/treetime/src/coalescent/integration.rs)
with a constant $T_c = 1$, rather than a dedicated integrator.

## Background

For constant $T_c$ the Kingman likelihood factorizes over intervals of constant
lineage count $k$, collapsing to

$$
L(T_c) \propto T_c^{-M}\,\exp(-I/T_c),
\qquad
I = \int \tfrac{k(k-1)}{2}\,dt,
$$

with $M$ the total number of merger events. $\partial_{T_c}\log L = 0$ yields the
MLE in closed form,

$$
T_c^{*} = I / M.
$$

$I$ and $M$ are accumulated per edge (see [features/coalescent.md](../features/coalescent.md)),
so they respect the same bad-branch exclusion as the likelihood and extend to the
skyline case.

## Rationale

- **Exact.** No bracket, tolerance, or iteration count. v0's bracket
  `[-20, 2]` in $\log T_c$ (roughly $[2\times10^{-9}, 7.4]$ linear) silently caps
  the estimate; the analytic form has no bounds artifact.
- **Cheaper.** One pass over the edges versus ~dozens of full-tree likelihood
  re-evaluations, each of which rebuilds the merger-rate integral.
- **Simpler to maintain.** Removes the `argmin` `BrentOpt` executor, the
  `TcCostFunction` cost-function plumbing, and the optimization observer from this
  path.
- **Result unchanged on real trees.** Within any edge's span $k \ge 2$, so the
  `max(0.5, k-1)` clamp that the likelihood path applies is never active there.
  Reusing `compute_integral_merger_rate` at $T_c = 1$ therefore yields exactly the
  textbook $\int (k-1)/2\,dt$ over edge endpoints and matches what the numerical
  optimizer converged to. The clamp in v0 was a guard for numerical integration,
  not a modeling choice.

## Options considered

1. **Keep Brent, port faithfully.** Rejected: retains iteration cost and bracket
   artifact for a problem with a one-line exact solution.
2. **Analytic, global integral of $k(k-1)/2$.** Correct for a clean tree but
   diverges from the per-edge likelihood once bad branches are excluded, and does
   not carry over to piecewise-constant $T_c$. Rejected in favor of the per-edge
   accumulation.
3. **Analytic, per-edge accumulation (chosen).** Matches the likelihood's edge
   set exactly and generalizes to the skyline.

## Impact

- Constant-$T_c$ estimates are numerically indistinguishable from v0 on
  well-behaved trees, and better-behaved on degenerate/extreme-timescale trees
  where Brent's bracket bit.
- The skyline path (piecewise-constant $T_c$) still uses numerical optimization
  (Nelder–Mead) and is unaffected.
- Fallback on a degenerate tree ($M \le 0$ or non-finite $I$) uses the initial
  $T_c$; the hardcoded initial value remains tracked in
  [issues/N-coalescent-initial-tc-hardcoded-fallback.md](../issues/N-coalescent-initial-tc-hardcoded-fallback.md).
