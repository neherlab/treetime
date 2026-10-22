# Coalescent

Kingman coalescent prior over node times, in decimal calendar coordinates. The
expected lineage counts `LineageExpectation` (see [algo/coalescent-contribution-refactor.md](../algo/coalescent-contribution-refactor.md))
supply the mean-field prior during timetree inference, the statistics of the $T_c$
fit, the polytomy merger rates, the reported likelihood and the convergence trace.

## Coalescent model

- [x] Calendar-coordinate expected lineage count $\bar k(t)=\sum_u w_u S_u(t)$ from node posteriors (or, in the init phase, Gaussian node marginals)
- [x] Mean-field prior per node, local in its own time
- [x] Per-lineage merger rate $\kappa(t)=\max(1/2,\bar k-1)/(2T_c)$ for polytomy resolution
- [x] Expected statistics per segment: pair-rate integral $I_i=\int \tfrac12(\bar k^2-\bar k+V)\,dt$ and mergers $M_i=\sum_v m_v P(t_v\in i)$, each node with its actual merger count
- [x] Bad branches and undated nodes count as surviving at all times
- [ ] v0 parity of the reported likelihood: v1 reports the expected complete-data log-likelihood, v0 the edge-based total likelihood; decision file pending (phase 5 of [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md))

## Time scale $T_c$

- [x] Constant $T_c$ (`--coalescent <TC>`): a bare number is in substitutions per site, the unit of v0 and augur, and is converted to years with the clock rate at the start of every refinement round and after the final clock estimate, as v0 divides by the clock rate ([packages/legacy/treetime/treetime/wrappers.py#L536-L538](../../packages/legacy/treetime/treetime/wrappers.py#L536-L538)); a number with the unit `y`, `w` or `d` is a time in years, weeks or days and stays fixed
- [x] Optimal constant $T_c$ (`--coalescent opt`) and skyline (`--coalescent skyline`): one setting with the fixed time scale, as in v0 and augur
- [x] **Optimal constant $T_c$ (analytic)** — see below
- [x] **Skyline (piecewise-constant $T_c$, analytic convex solve)** — see below
- [x] No hardcoded initial/fallback $T_c$ — the analytic solve needs no seed and
  cannot fail numerically; it errors only on a tree that is degenerate for the
  coalescent (no time span or no mergers), which the pipeline propagates so the run
  stops with a clear message rather than substituting an invented timescale

## Optimal constant $T_c$

A constant $T_c$ is exactly the one-segment skyline, so
`estimate_coalescent_tc()` in
[`packages/treetime/src/timetree/coalescent_timescale.rs#L40`](../../packages/treetime/src/timetree/coalescent_timescale.rs#L40)
calls `optimize_skyline` with `n_points = 1` for the constant mode and reports its single
segment's $T_c$ and likelihood. The derivation below explains why that one-segment
solve is the closed form $T_c^{*}=I/M$; the I/M accumulation, likelihood reporting,
and degenerate-tree error are the shared skyline machinery (see
[Skyline](#skyline-piecewise-constant-t_c)).

The constant-$T_c$ Kingman likelihood factorizes over intervals of constant
lineage count $k$. Each merger event contributes a factor $k(k-1)/(2T_c)$ and each
interval of duration $\Delta t$ a survival factor $\exp\!\big(-\Delta t\,k(k-1)/(2T_c)\big)$.
Taking the product over the whole tree collapses to a closed form in $T_c$:

$$
L(T_c) \propto T_c^{-M}\,\exp\!\left(-\frac{I}{T_c}\right),
\qquad
\log L = -M\ln T_c - \frac{I}{T_c} + \text{const},
$$

where

$$
I = \int \frac{k(k-1)}{2}\,dt
\qquad\text{(pairwise-rate integral, } T_c \text{ factored out),}
$$

and $M$ is the total number of merger events. Setting $\partial_{T_c}\log L = 0$
gives the maximum-likelihood estimate directly:

$$
\boxed{\,T_c^{*} = I / M\,}.
$$

No iterative search is needed. This replaces v0's Brent minimization of
`-total_LH()` over $\log T_c$; see [decisions/coalescent-analytic-tc-optimization.md](../decisions/coalescent-analytic-tc-optimization.md).

### Per-segment statistics (shared with the skyline)

The sufficient statistics are the expected $I_i$ and $M_i$ of the node posteriors
(`LineageExpectation::segment_statistics`); the constant case is the single segment.
With $E[k(k-1)] = \bar k^2 - \bar k + V$ for independent node times, $I_i$ reduces to
the textbook $\int k(k-1)/2\,dt$ when the node times are points, and $M$ is the total
merger count, $N-1$ for a clean bifurcating tree. The integral runs on the expectation's
grid of 4001 points, so a step of $k$ at point times costs half a grid cell. The segments
span from where $\bar k$ rises half a lineage above its value before every node, near the
root, to the latest tip (`LineageExpectation::span`).

### Reported likelihood and degenerate-tree error

The reported coalescent log-likelihood at the optimum is the expected complete-data
log-likelihood $-\sum_i I_i/T_{c,i} - \sum_i M_i \ln T_{c,i}$, the objective of the solve
(`fn expected_log_likelihood`); the convergence trace evaluates the same quantity for the
posteriors of each round. The closed form $T_c^{*}=I/M$ has no numerical failure mode; the
only way it errors is a tree that is degenerate for the coalescent — $M \le 0$ (no
mergers) or $I \le 0$ / non-finite (no time span). In that case the shared skyline path
returns an error (no fallback timescale), which the pipeline propagates to stop the run.

## Skyline (piecewise-constant $T_c$)

Implemented in [`packages/treetime/src/coalescent/skyline.rs`](../../packages/treetime/src/coalescent/skyline.rs).

The time span is split into `n_points` equal-width segments; $T_c$ is constant
within each. Writing $z_i = \ln T_{c,i}$ (so the coalescence rate is
$1/T_{c,i} = e^{-z_i}$), the negative log-likelihood plus a smoothness penalty is

$$
C(z) = \sum_i \big(I_i\,e^{-z_i} + M_i z_i\big) + \frac{\gamma}{2}\sum_i (z_{i+1}-z_i)^2,
$$

with per-segment $I_i = \int_{\text{seg }i} k(k-1)/2\,dt$, merger count $M_i$, and
stiffness $\gamma$. Modeling $z = \ln T_c$ makes the penalty **scale-independent** —
it charges squared *log-fold-changes* $z_{i+1}-z_i = \ln(T_{c,i+1}/T_{c,i})$ — and
guarantees $T_c = e^z > 0$ with no constraint. Every term is convex in $z$, so $C$
has a **unique minimizer** (a Poisson-likelihood / Gaussian-smoothing MAP problem —
the frequentist analog of the Bayesian skygrid). It is found by **Newton's method on
the symmetric tridiagonal Hessian** ($O(n)$ Thomas solves), warm-started from the
decoupled per-segment optimum $z_i = \ln(I_i/M_i)$ and globalized with an Armijo
backtracking line search.

- [x] Change of variable to $z = \ln T_c$ (convex objective, unique optimum,
  scale-independent smoothing, positivity for free)
- [x] Expected per-segment $I_i$, $M_i$ of the node posteriors; the reported LH is the
  objective the solve maximizes
- [x] Newton solve with Armijo line search (no positivity constraint needed)
- [x] Curvature-based confidence bands from the full inverse Hessian in $\ln T_c$ coordinates, including adjacent-segment covariance from the stiffness penalty ([kb/decisions/coalescent-skyline-hessian-confidence-bands.md](../decisions/coalescent-skyline-hessian-confidence-bands.md))
- [x] **Equal-width segment boundaries** — uniform spacing gives the stiffness a
  clean, grid-independent meaning and is robust to tied merger times. Merger-sparse
  regions may leave empty segments ($M_i=0$) whose $z_i$ is pinned by the smoothing
  prior, so `stiffness > 0` is required (and enforced) for `n_points > 1` to keep the
  Hessian positive-definite and every $z_i$ finite (see
  [decisions/coalescent-skyline-convex-log-tc.md](../decisions/coalescent-skyline-convex-log-tc.md))
- [x] Piecewise-constant $T_c(t)$ output; times outside the grid clamp to the
  first/last segment

This replaces the previous Nelder–Mead search over piecewise-linear $\log T_c$
knots; see [decisions/coalescent-skyline-convex-log-tc.md](../decisions/coalescent-skyline-convex-log-tc.md).
