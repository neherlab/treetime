# Skyline Tc solved analytically in 1/Tc instead of Nelder-Mead over log Tc

v1 estimates the piecewise-constant skyline $T_c(t)$ by an exact convex
optimization in the inverse time scale $x_i = 1/T_{c,i}$. The previous v1
implementation (and v0) used a derivative-free / general-purpose numerical search
over $\log T_c$ knots.

**Type**: Implementation change (new parametrization and solver; also diverges
from v0's optimizer, extending the constant-$T_c$ decision to the skyline).

**Previous v1**: `optimize_skyline()` used `argmin`'s `NelderMead` over a
piecewise-_linear_ $\log T_c$ grid with a stiffness penalty on adjacent $\log T_c$
and a boundary penalty, plus a hand-built initial simplex.

**v0**: `MergerModel.optimize_skyline()` at
[`packages/legacy/treetime/treetime/merger_models.py#L281-L318`](../../packages/legacy/treetime/treetime/merger_models.py) uses
`scipy.optimize.minimize(method='SLSQP')`.

**v1 now**: `optimize_skyline()` at
[`packages/treetime/src/coalescent/skyline.rs`](../../packages/treetime/src/coalescent/skyline.rs).

## Background

Splitting the tree into segments with constant $T_c$, the constant-$T_c$ likelihood
(see [coalescent-analytic-tc-optimization.md](coalescent-analytic-tc-optimization.md))
applies per segment. Writing $x_i = 1/T_{c,i}$ and adding a smoothness penalty
$\frac{\gamma}{2}\sum_i (x_{i+1}-x_i)^2$, the cost to minimize is

$$
C(x) = \sum_i \big(I_i x_i - M_i \ln x_i\big) + \frac{\gamma}{2}\sum_i (x_{i+1}-x_i)^2,
$$

with $I_i = \int_{\text{seg }i} k(k-1)/2\,dt$ and merger count $M_i$. Each term is
convex on $x_i>0$, so $C$ has a unique positive minimizer, found by Newton's method
on the symmetric tridiagonal Hessian (Thomas algorithm, $O(n)$ per step),
warm-started from the decoupled optimum $x_i = M_i/I_i$ and damped to keep
$x_i>0$. This is a Poisson-likelihood / Gaussian-smoothing MAP estimate — the
frequentist analog of the Bayesian GMRF skygrid.

## Rationale

- **Convex with a unique optimum.** In $x = 1/T_c$ the log-likelihood is concave
  and the penalty is a PSD quadratic, so there are no local optima and no simplex
  or bracket to tune. Optimizing $T_c$ (or $\log T_c$) directly, as before, is not
  convex.
- **Exact and cheap.** A few $O(n)$ Newton steps versus many derivative-free
  objective evaluations, each rebuilding the merger-rate integral.
- **Self-consistent likelihood.** $I_i$ and $M_i$ use the same interval-midpoint
  and node-time attribution as `CoalescentModel`, so the analytic optimum is the
  true maximizer of the model-evaluated likelihood, and the reported LH matches
  `compute_coalescent_total_lh`.
- **Simpler.** Removes the Nelder-Mead executor, the initial-simplex construction,
  the boundary penalty, and the piecewise-linear knot interpolation.

## Merger-quantile segment boundaries

Segment boundaries are the quantiles of the merger times (outer boundaries pinned
to the tree's time span). This is deliberate: with equal-time boundaries the
boundary segment between the youngest coalescence and the tips has lineages
($I_i>0$) but no mergers ($M_i=0$). With $M_i=0$ the $-M_i\ln x_i$ barrier
vanishes and the segment's optimum sits at $x_i\to 0$, i.e. $T_c\to\infty$ (only
weakly restrained by smoothing). Quantile boundaries let every segment own mergers
so $M_i>0$ keeps every $x_i$ bounded away from zero. On ebola/20 this changed a
divergent most-recent segment ($T_c\approx 3\times10^{150}$) into a finite, smoothly
varying trajectory.

This holds only when the requested segment count does not exceed the number of
distinct merger times. When it does, the clamped quantile construction emits
duplicate (zero-width) or zero-exposure segments, for which $I_i x - M_i \ln x$ is
unbounded below. `optimize_skyline` now rejects such an over-segmented grid up front
(and validates strictly increasing boundaries and positive per-segment exposure as
backstops) rather than masking the degeneracy with the pooled warm-start fallback,
so the guarantee is enforced instead of assumed.

## Options considered

1. **Port v0's SLSQP over log Tc.** Rejected: retains non-convex parametrization
   and optimizer-dependent local optima; the constant-$T_c$ case already showed the
   numerical search landing on wrong values.
2. **Newton in $x=1/T_c$, equal-time boundaries.** Convex and exact, but empty
   boundary segments drive $T_c\to\infty$.
3. **Newton in $x=1/T_c$, merger-quantile boundaries (chosen).** Convex, exact, and
   structurally free of empty-segment blow-up.

## Impact

- Skyline estimates no longer match v0 numerically; they are the exact optimum of
  the (regularized) piecewise-constant Kingman likelihood.
- The stiffness parameter now penalizes adjacent $1/T_c$ differences rather than
  $\log T_c$ differences, so its scale differs from the previous implementation.
- Supersedes the Nelder-Mead optimizer and simplex-initialization concerns in
  [issues/M-timetree-skyline-nelder-mead-optimizer.md](../issues/M-timetree-skyline-nelder-mead-optimizer.md)
  and [issues/N-coalescent-skyline-simplex-initialization-undecided.md](../issues/N-coalescent-skyline-simplex-initialization-undecided.md).
  Boundary out-of-domain evaluation still uses constant extrapolation
  ([issues/N-coalescent-skyline-extrapolation-policy-undecided.md](../issues/N-coalescent-skyline-extrapolation-policy-undecided.md)).
