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
- [x] Skyline (piecewise-constant $T_c$, Nelder–Mead over log-$T_c$ knots)
- [/] Data-derived initial/fallback $T_c$ (hardcoded constant, see [issues/N-coalescent-initial-tc-hardcoded-fallback.md](../issues/N-coalescent-initial-tc-hardcoded-fallback.md))

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
- $M = \sum_{\text{edges}} (n_\text{children}-1)/n_\text{children}$. Summed over a
  node's $n_\text{children}$ edges this is the node's merger count
  $n_\text{children}-1$; over a clean bifurcating tree, $M = N-1$.

`compute_integral_merger_rate` applies the `max(0.5, k-1)` clamp used on the
likelihood path, but within any edge's span the lineage count satisfies $k \ge 2$
(the clamp only activates where $k<1.5$, i.e. above the root, where no edge
exists). So evaluating it at $T_c = 1$ over edge endpoints yields exactly the
textbook $\int (k-1)/2\,dt$, and no separate unclamped integrator is needed.

### Reported likelihood and fallback

`optimize_tc` still reports the total coalescent log-likelihood at $T_c^{*}$ by
evaluating the shared `CoalescentModel`, so the value matches
`compute_coalescent_total_lh` bit-for-bit. If $M \le 0$ or $I$/$T_c^{*}$ are
non-finite (degenerate tree), it falls back to the caller-supplied initial $T_c$
and reports `success = false`.
