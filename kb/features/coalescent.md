# Coalescent

Kingman coalescent prior over node times, in decimal calendar coordinates. The
shared `CoalescentModel` (see [algo/coalescent-contribution-refactor.md](../algo/coalescent-contribution-refactor.md))
supplies per-node contributions during timetree inference and endpoint-derived
per-edge contributions for the whole-tree likelihood and for $T_c$ optimization.

## Coalescent model

- [x] Calendar-coordinate lineage count $k(t)$ (piecewise constant)
- [x] Per-lineage merger rate $\kappa(t)=\max(0.5,k-1)/(2T_c)$
- [x] Total pairwise merger rate $\lambda(t)=k(k-1)/(2T_c)$
- [x] Expected-merger integral $H(t)=\int_t^P \kappa(s)\,ds$
- [x] Per-edge contribution with actual parent multiplicity (see [decisions/coalescent-total-lh-actual-multiplicity.md](../decisions/coalescent-total-lh-actual-multiplicity.md))
- [x] Bad branches excluded from edge collection

## Time scale $T_c$

- [x] Constant $T_c$
- [x] **Optimal constant $T_c$ (analytic)** — see below
- [x] **Skyline (piecewise-constant $T_c$, analytic convex solve)** — see below
- [x] No hardcoded initial/fallback $T_c$ — the analytic solve needs no seed and
  cannot fail numerically; it errors only on a tree that is degenerate for the
  coalescent (no time span or no mergers), which the pipeline propagates so the run
  stops with a clear message rather than substituting an invented timescale (see
  [issues/N-coalescent-initial-tc-hardcoded-fallback.md](../issues/N-coalescent-initial-tc-hardcoded-fallback.md))

## Optimal constant $T_c$

Implemented in [`packages/treetime/src/coalescent/optimize_tc.rs`](../../packages/treetime/src/coalescent/optimize_tc.rs).

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

### Per-edge accumulation

Both $I$ and $M$ are accumulated per edge rather than as global quantities, so
they stay consistent with the bad branches that `collect_coalescent_edges`
excludes, and so the expressions generalize to piecewise-constant $T_c$ (where
per-interval merger counts matter):

- $I = \sum_{\text{edges}} \big[H_0(t_\text{parent}) - H_0(t_\text{child})\big]$,
  where $H_0(t)=\int_t^P (k-1)/2\,dt$ is the $T_c$-independent per-lineage
  integral. It is obtained by reusing
  [`compute_integral_merger_rate`](../../packages/treetime/src/coalescent/integration.rs)
  with a constant $T_c = 1$ rather than a dedicated function. Summed over the $k$
  edges spanning an interval, the per-lineage integrand $(k-1)/2$ reproduces the
  pairwise integrand $k(k-1)/2$.
- $M = \sum_{\text{edges}} (n_\text{siblings}-1)/n_\text{siblings}$, where
  $n_\text{siblings}$ is the number of children of the edge's parent node (this
  child plus its siblings). Summed over a node's $n_\text{siblings}$ edges this is
  the node's merger count $n_\text{siblings}-1$; over a clean bifurcating tree,
  $M = N-1$.

`compute_integral_merger_rate` applies the `max(0.5, k-1)` clamp used on the
likelihood path, but within any edge's span the lineage count satisfies $k \ge 2$
(the clamp only activates where $k<1.5$, i.e. above the root, where no edge
exists). So evaluating it at $T_c = 1$ over edge endpoints yields exactly the
textbook $\int (k-1)/2\,dt$, and no separate unclamped integrator is needed.

### Reported likelihood and degenerate-tree error

`optimize_tc` reports the total coalescent log-likelihood at $T_c^{*}$ by
evaluating the shared `CoalescentModel`, so the value matches
`compute_coalescent_total_lh` bit-for-bit. The closed form $T_c^{*}=I/M$ has no
numerical failure mode; the only way it errors is a tree that is degenerate for the
coalescent — $M \le 0$ (no mergers) or $I \le 0$ / non-finite (no time span). In
that case it returns an error (no fallback timescale), which the pipeline
propagates to stop the run.

## Skyline (piecewise-constant $T_c$)

Implemented in [`packages/treetime/src/coalescent/skyline.rs`](../../packages/treetime/src/coalescent/skyline.rs).

The time span is split into `n_points` segments; $T_c$ is constant within each.
Writing $z_i = \ln T_{c,i}$ (so the coalescence rate is $1/T_{c,i} = e^{-z_i}$), the
negative log-likelihood plus a smoothness penalty is

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
- [x] Analytic per-segment $I_i$, $M_i$ (same interval-midpoint / node-time
  conventions as `CoalescentModel`, so the optimum maximizes the model-evaluated
  likelihood and the reported LH matches `compute_coalescent_total_lh`)
- [x] Newton solve with Armijo line search (no positivity constraint needed)
- [x] **Merger-quantile segment boundaries** — boundaries fall at quantiles of the
  merger times so every segment (including the root and tip boundary segments)
  owns mergers; $M_i>0$ keeps the linear $M_i z_i$ term active, bounding $z_i$ from
  above and preventing empty boundary segments from collapsing to $T_c\to\infty$
  (see [decisions/coalescent-skyline-convex-log-tc.md](../decisions/coalescent-skyline-convex-log-tc.md))
- [x] Piecewise-constant $T_c(t)$ output; times outside the grid clamp to the
  first/last segment

This replaces the previous Nelder–Mead search over piecewise-linear $\log T_c$
knots; see [decisions/coalescent-skyline-convex-log-tc.md](../decisions/coalescent-skyline-convex-log-tc.md).
