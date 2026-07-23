# Skyline Tc solved by exact convex Newton in log Tc

v1 estimates the piecewise-constant skyline $T_c(t)$ by an exact convex
optimization in the log time scale $z_i = \ln T_{c,i}$. This replaced two earlier
approaches: the original v1 derivative-free search over $\log T_c$ knots, and an
intermediate exact convex solve in the inverse time scale $x_i = 1/T_{c,i}$.

**Type**: Implementation change (parametrization and solver; also diverges from
v0's optimizer, extending the constant-$T_c$ decision to the skyline).

**Original v1**: `optimize_skyline()` used `argmin`'s `NelderMead` over a
piecewise-*linear* $\log T_c$ grid with a stiffness penalty on adjacent $\log T_c$
and a boundary penalty, plus a hand-built initial simplex.

**Intermediate v1**: exact convex Newton in $x_i = 1/T_{c,i}$ minimizing
$\sum_i (I_i x_i - M_i \ln x_i) + \frac{\gamma}{2}\sum_i (x_{i+1}-x_i)^2$, with
positivity-preserving damped steps.

**v0**: `MergerModel.optimize_skyline()` at
[`packages/legacy/treetime/treetime/merger_models.py#L281-L318`](../../packages/legacy/treetime/treetime/merger_models.py) uses
`scipy.optimize.minimize(method='SLSQP')`.

**v1 now**: `optimize_skyline()` at
[`packages/treetime/src/coalescent/skyline.rs`](../../packages/treetime/src/coalescent/skyline.rs).

## Background

Splitting the tree into segments with constant $T_c$, the constant-$T_c$ likelihood
(see [coalescent-analytic-tc-optimization.md](coalescent-analytic-tc-optimization.md))
applies per segment. Writing $z_i = \ln T_{c,i}$ (so the coalescent rate is
$1/T_{c,i} = e^{-z_i}$) and penalizing squared log-fold-changes of $T_c$, the cost
to minimize is

$$
C(z) = \sum_i \big(I_i\,e^{-z_i} + M_i z_i\big) + \frac{\gamma}{2}\sum_i (z_{i+1}-z_i)^2,
$$

with $I_i = \int_{\text{seg }i} k(k-1)/2\,dt$ and merger count $M_i$. Each data term
is convex ($e^{-z_i}$ convex, $M_i z_i$ linear) and the penalty is a PSD quadratic,
so $C$ has a unique minimizer, found by Newton's method on the symmetric tridiagonal
Hessian (Thomas algorithm, $O(n)$ per step), warm-started from the decoupled optimum
$z_i = \ln(I_i/M_i)$ and globalized with an Armijo line search. This is a
Poisson-likelihood / Gaussian-smoothing MAP estimate — the frequentist analog of the
Bayesian GMRF skygrid, which is likewise formulated on $\log$ effective population
size.

## Rationale

- **Convex with a unique optimum.** In $z = \ln T_c$ the data term is convex and the
  penalty is a PSD quadratic, so there are no local optima and no simplex or bracket
  to tune. Optimizing $T_c$ directly is not convex.
- **Scale-independent smoothing.** The penalty charges $(z_{i+1}-z_i)^2 =
  \ln(T_{c,i+1}/T_{c,i})^2$, i.e. squared log-fold-changes, so the stiffness is
  dimensionless and its effect does not depend on the absolute scale of $T_c$ (or the
  time unit). In the earlier $x = 1/T_c$ parametrization the penalty
  $(x_{i+1}-x_i)^2$ carried units of time$^{-2}$, so the same stiffness meant
  different smoothing at different scales.
- **Positivity for free.** $T_c = e^z > 0$ for any real $z$, so the Newton step needs
  no positivity cap — only an Armijo line search for global convergence (the data
  term is no longer quadratic in the optimization variable, so Newton is iterative).
- **Exact and cheap.** A few $O(n)$ Newton steps versus many derivative-free
  objective evaluations, each rebuilding the merger-rate integral.
- **Self-consistent likelihood.** $I_i$ and $M_i$ use the same interval-midpoint and
  node-time attribution as `CoalescentModel`, so the analytic optimum is the true
  maximizer of the model-evaluated likelihood, and the reported LH matches
  `compute_coalescent_total_lh`.
- **Simpler.** Removes the Nelder-Mead executor, the initial-simplex construction,
  the boundary penalty, and the piecewise-linear knot interpolation.

## Merger-quantile segment boundaries

Segment boundaries are the quantiles of the merger times (outer boundaries pinned
to the tree's time span). This is deliberate: with equal-time boundaries the
boundary segment between the youngest coalescence and the tips has lineages
($I_i>0$) but no mergers ($M_i=0$). With $M_i=0$ the data term reduces to
$I_i e^{-z_i}$, whose gradient $-I_i e^{-z_i}$ is negative everywhere, so the
segment's optimum runs to $z_i\to+\infty$, i.e. $T_c\to\infty$ (only weakly
restrained by smoothing). Quantile boundaries guarantee every segment — including
both boundary segments — owns mergers, so $M_i>0$ keeps the linear $M_i z_i$ term
active and bounds $z_i$ from above. (In the earlier $x = 1/T_c$ parametrization the
analogous guard was the $-M_i\ln x_i$ log-barrier keeping $x_i$ away from zero.) On
ebola/20 quantile boundaries changed a divergent most-recent segment
($T_c\approx 3\times10^{150}$) into a finite, smoothly varying trajectory.

## Options considered

1. **Port v0's SLSQP over log Tc.** Rejected: retains a general-purpose numerical
   search with optimizer-dependent local optima; the constant-$T_c$ case already
   showed the numerical search landing on wrong values.
2. **Newton in $x=1/T_c$ (intermediate v1).** Convex and exact with a tridiagonal
   Hessian, but the smoothing penalty $(x_{i+1}-x_i)^2$ is scale-dependent and the
   step needs a positivity cap to keep $x_i>0$.
3. **Newton in $z=\ln T_c$, merger-quantile boundaries (chosen).** Convex and exact
   with the same tridiagonal Hessian, plus scale-independent smoothing and
   unconstrained positivity. The data term becomes non-quadratic (exponential), so
   Newton iterates under an Armijo line search rather than converging in one step.

## Impact

- Skyline estimates no longer match v0 numerically; they are the exact optimum of
  the (regularized) piecewise-constant Kingman likelihood.
- The stiffness parameter penalizes adjacent $\ln T_c$ differences
  (log-fold-changes), so it is dimensionless and scale-independent. Its numeric
  scale differs from the intermediate $1/T_c$ implementation; the CLI default is
  retained at `2.0` and reinterpreted on the log scale.
- Supersedes the Nelder-Mead optimizer and simplex-initialization concerns in
  [issues/M-timetree-skyline-nelder-mead-optimizer.md](../issues/M-timetree-skyline-nelder-mead-optimizer.md)
  and [issues/N-coalescent-skyline-simplex-initialization-undecided.md](../issues/N-coalescent-skyline-simplex-initialization-undecided.md).
  Boundary out-of-domain evaluation still uses constant extrapolation
  ([issues/N-coalescent-skyline-extrapolation-policy-undecided.md](../issues/N-coalescent-skyline-extrapolation-policy-undecided.md)).
