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
  - $m$: fitted linear coefficient (field `slope`); $m = 0$ when $b > 0$, $m = s$ when $b = 0$

  Fitted from the innermost grid points [`packages/treetime-grid/src/hard_approach_law.rs`](../../packages/treetime-grid/src/hard_approach_law.rs). Same functional form as truth.

- **Part C** ([`kb/proposals/distribution-log-space-and-hard-soft-boundaries.md#L158-L169`](../proposals/distribution-log-space-and-hard-soft-boundaries.md#L158-L169)):

  $y = a - b\ln\delta$

  No linear term. Its $b = 0$ case is flat ($y = a$), yet the proposal's own truth states $y = \mu t + C$ (linear) for $n = 0$ ([`kb/proposals/distribution-log-space-and-hard-soft-boundaries.md#L149`](../proposals/distribution-log-space-and-hard-soft-boundaries.md#L149), [`kb/proposals/distribution-log-space-and-hard-soft-boundaries.md#L156`](../proposals/distribution-log-space-and-hard-soft-boundaries.md#L156)). Internally inconsistent.

## Choice 1: drop the linear term for divergent branches ($b > 0$)

- The fit sets $m = 0$ [`packages/treetime-grid/src/hard_approach_law.rs#L105-L121`](../../packages/treetime-grid/src/hard_approach_law.rs#L105-L121), omitting $s\,\delta$.
- Bounded: over the sub-grid gap $\delta < 0.01\,\ell$ (with $\ell$ one mutation's branch length), so $|s\,\delta| < 0.01$ in neg-log.
- Fitting a slope next to a $\ln\delta$ divergence from five points is numerically unstable.
- The $b = 0$ case keeps $m = s$ (linear), reproducing the mode on the boundary.

**Open question**: accept the dropped term as a bounded approximation, or seed $m = s$ analytically at construction (the rate is known there), making the law exact without a fit.

## Choice 2: fold indel counts into one exponent

- The fit reads stored ordinates that already include the indel term, so it recovers a single $b \approx n + k$ and cannot separate substitution from indel events.
- Consequence of adding the Poisson indel term to the branch-length likelihood, a v1 decision. The boundary divergence steepens with indel events exactly as with substitutions.

**Open question**: confirm that a single combined exponent is the intended treatment, or that substitution and indel events should enter the boundary law differently.

## Impact

- No demonstrated wrong result.
- The approach law governs only $[t_\text{hard}, t_\text{first})$, whose mass is small. The grid interior evaluates the exact likelihood.
- Open science questions. No `kb/decisions` entry until reviewed.

## Related

- [M-grid-hard-approach-fit-misclassifies-falling-linear.md](M-grid-hard-approach-fit-misclassifies-falling-linear.md): separate fit-robustness concern for the same law.
- [kb/proposals/distribution-log-space-and-hard-soft-boundaries.md#L140](../proposals/distribution-log-space-and-hard-soft-boundaries.md#L140): Part C, which omits the linear term in both regimes.
