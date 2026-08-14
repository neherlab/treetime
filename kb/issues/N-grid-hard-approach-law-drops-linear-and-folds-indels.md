# Hard-approach boundary law drops the linear term for divergent branches and folds indel counts into the exponent

The branch-length hard-approach law makes two modeling choices near the `t = 0` boundary that need review by the science team. Neither has a demonstrated wrong result; both are bounded approximations, but they are scientific choices, not mechanical ones.

## Background

The branch-length neg-log likelihood near the boundary is, to leading order, the Gamma kernel `y(t) = -(n + k)·ln t + (mu + r)·t + C`, where `n` is the substitution mutation count, `k` the indel count, `mu` the aggregate substitution rate, and `r` the indel rate. This is confirmed analytically (transition-probability product plus a Poisson indel term) and matches the stored grid, which sums the substitution log-likelihood and the Poisson indel term [`packages/treetime/src/optimize/likelihood.rs#L89-L98`](../../packages/treetime/src/optimize/likelihood.rs#L89-L98).

`HardApproachLaw` fits `y(t) = a - b·ln|t - t_hard| + slope·(t - t_hard)` from the innermost grid points [`packages/treetime-grid/src/hard_approach_law.rs#L65-L139`](../../packages/treetime-grid/src/hard_approach_law.rs#L65-L139).

## Choice 1: drop the linear term for divergent branches

For a branch with events (`b > 0`) the fit sets `slope = 0` [`packages/treetime-grid/src/hard_approach_law.rs#L105-L121`](../../packages/treetime-grid/src/hard_approach_law.rs#L105-L121), so the law omits the `(mu + r)·t` term. Over the sub-grid gap the omitted term is small (`t < 0.01·one_mutation`, so the term is under about `0.01` in neg-log), and fitting a slope next to a logarithmic divergence from five points is numerically unstable. For a zero-event branch (`b = 0`) the linear term is kept exactly, which reproduces the mode on the boundary.

Open question: accept the dropped term as a bounded approximation, or seed `slope = mu + r` analytically at construction, where the rate is known, so the law is exact without a fit.

## Choice 2: fold indel counts into one boundary exponent

The fit reads stored ordinates that already include the indel term, so it recovers a single combined exponent `b ≈ n + k` and cannot separate substitution from indel events. This follows from adding the Poisson indel term to the branch-length likelihood, itself a v1 decision. The boundary divergence therefore steepens with indel events exactly as with substitutions.

Open question: confirm that a single combined power-law exponent at the boundary is the intended scientific treatment of substitution and indel events, or that they should enter the boundary law differently.

## Impact

No demonstrated wrong result. The approach law governs only the region between the boundary and the first grid point, whose mass is small; the grid interior evaluates the exact likelihood. These are open scientific questions to settle before recording a `kb/decisions` entry for the hard-approach law.

## Related

- [M-grid-hard-approach-fit-misclassifies-falling-linear.md](M-grid-hard-approach-fit-misclassifies-falling-linear.md): a separate fit-robustness concern for the same law.
- [kb/proposals/distribution-log-space-and-hard-soft-boundaries.md](../proposals/distribution-log-space-and-hard-soft-boundaries.md): Part C, which specifies only `a - b·log(t - t_hard)` and omits the linear term in both regimes.
