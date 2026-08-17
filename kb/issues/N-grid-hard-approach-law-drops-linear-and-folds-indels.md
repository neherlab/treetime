# Hard-approach boundary law drops the linear term for divergent branches and folds indel counts into the exponent

`struct HardApproachLaw` encodes two modeling choices near the $t = 0$ boundary that need science review.

## The three laws

- **truth** (leading order):

  $y = -b\ln\delta + s\,\delta + C$

  where:
  - $\delta = t - t_\text{hard}$: distance from the hard boundary ($t_\text{hard} = 0$ for branch length)
  - $b = n + k$: power-law exponent; substitution count $n$ plus indel count $k$
  - $s = \mu + r$: linear coefficient; substitution rate $\mu$ plus indel rate $r$
  - $C$: normalization constant

  Confirmed analytically (transition-probability product plus a Poisson indel term) and against the stored grid, which sums the substitution and Poisson indel log-likelihoods [`packages/treetime/src/optimize/likelihood.rs#L89-L98`](../../packages/treetime/src/optimize/likelihood.rs#L89-L98).

- **code** (`struct HardApproachLaw`):

  $y = y_\text{edge} - b\ln(\delta/\delta_\text{edge}) + m\,\delta$

  where:
  - $y_\text{edge}$: the live neg-log grid edge ordinate, read on evaluation (edge-relative; no stored anchor, only the fixed boundary location $t_\text{hard}$)
  - $m$: linear coefficient (field `slope`); $m = 0$ when $b > 0$ (divergent), $m = (y_\text{edge} - y_\text{hard}) / (t_\text{edge} - t_\text{hard})$ when $b = 0$ (finite)

  The regime is set by the evaluated boundary ordinate $y_\text{hard}$: infinite is divergent (fit $b$ from the innermost grid points, $m = 0$); finite is the boundary-anchored line ($b = 0$, exact $m$ from the boundary to the grid edge) [`packages/treetime-grid/src/hard_approach_law.rs`](../../packages/treetime-grid/src/hard_approach_law.rs). Same functional form as truth.

- **Part C** ([`kb/proposals/distribution-log-space-and-hard-soft-boundaries.md#L158-L169`](../proposals/distribution-log-space-and-hard-soft-boundaries.md#L158-L169)):

  $y = a - b\ln\delta$

  No linear term. Its $b = 0$ case is flat ($y = a$), yet the proposal's own truth states $y = \mu t + C$ (linear) for $n = 0$ ([`kb/proposals/distribution-log-space-and-hard-soft-boundaries.md#L149`](../proposals/distribution-log-space-and-hard-soft-boundaries.md#L149), [`kb/proposals/distribution-log-space-and-hard-soft-boundaries.md#L156`](../proposals/distribution-log-space-and-hard-soft-boundaries.md#L156)). Internally inconsistent.

## Choice 1: drop the linear term for divergent branches ($b > 0$)

- The divergent branch sets $m = 0$ [`packages/treetime-grid/src/hard_approach_law.rs`](../../packages/treetime-grid/src/hard_approach_law.rs), omitting $s\,\delta$.
- Bounded: over the sub-grid gap $\delta < 0.01\,\ell$ (with $\ell$ one mutation's branch length), so $|s\,\delta| < 0.01$ in neg-log.
- The finite branch ($b = 0$) is boundary-anchored: $m$ is the exact line from the evaluated boundary ordinate to the grid edge, reproducing the mode on the boundary without a fit.

**Open question**: for the divergent branch, accept the dropped linear term as a bounded approximation, or seed $m = s$ analytically at construction (the rate is known there), making the law exact.

## Choice 2: fold indel counts into one exponent

- The fit reads stored ordinates that already include the indel term, so it recovers a single $b \approx n + k$ and cannot separate substitution from indel events.
- Consequence of adding the Poisson indel term to the branch-length likelihood, a v1 decision. The boundary divergence steepens with indel events exactly as with substitutions.

**Open question**: confirm that a single combined exponent is the intended treatment, or that substitution and indel events should enter the boundary law differently.

## Impact

- No demonstrated wrong result.
- The approach law governs only $[t_\text{hard}, t_\text{first})$, whose mass is small. The grid interior evaluates the exact likelihood.
- Open science questions. No `kb/decisions` entry until reviewed.

## Related

- [kb/proposals/distribution-log-space-and-hard-soft-boundaries.md#L140](../proposals/distribution-log-space-and-hard-soft-boundaries.md#L140): Part C, which omits the linear term in both regimes.
