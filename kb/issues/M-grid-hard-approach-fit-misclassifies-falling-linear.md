# Hard-approach fit misclassifies a falling linear neg-log density as a power law

`fn HardApproachLaw::fit()` on a strictly decreasing linear neg-log grid returns a power-law exponent `b > 0` instead of the linear case `b = 0`, so the linear refit that would recover the `slope` never runs and the boundary shape is misrepresented.

## Symptom and reproduction

`test_hard_approach_law_fit_recovers_linear::case_3_falling` (`slope = -2.0`, `intercept = 4.0`) fails: it expects `b = 0` but the fit returns `b ≈ 0.516`. The stored grid is the exact line `y(t) = -2 t + 4` over `t ∈ [0.1, 1.0]`, so the correct recovery is `b = 0`, `slope = -2.0`, `a = 4.0`. The rising, flat, and steep-rising cases all pass; only the falling (negative-slope) case fails.

## Root cause

`fit` classifies the boundary shape in stage 1 by least-squares regression of `(ln|t - t_hard|, y)`, taking the slope as `-b` [`packages/treetime-grid/src/hard_approach_law.rs#L84-L101`](../../packages/treetime-grid/src/hard_approach_law.rs#L84-L101). Over a fit window bounded away from `t_hard`, a falling line `y` decreases as `ln|t - t_hard|` increases, so the regression slope is negative and `b` comes out spuriously positive (here `≈ 0.516`). Because the linear refit only runs when `b <= B_THRESHOLD` (`0.01`) [`packages/treetime-grid/src/hard_approach_law.rs#L103-L128`](../../packages/treetime-grid/src/hard_approach_law.rs#L103-L128), a genuinely linear but falling density never reaches it and is reported as a power law. A rising line gives a negative raw `b` that clamps to `0` and passes; a falling line does not.

## Impact

Latent today: no production path fits a `HardApproachLaw` yet, so no current output is wrong. It matters for the log-space boundary work, where the `n = 0` zero-mutation branch has a neg-log density that is finite and can fall toward the boundary; misclassifying it as a vanishing power law would misplace the boundary approach and any mode sitting on it.

## Fix approach

The `b >= 0` clamp cannot distinguish a falling line (spurious positive `b`) from a true power law. Detect linearity independently of the log-distance fit, for example by comparing the residuals of a linear-in-`t` fit against the linear-in-`ln|dt|` fit and selecting the better model, or by restricting the power-law branch to windows where `y` increases toward `t_hard`.
