# Hard-approach fit misclassifies a falling linear neg-log density as a power law

`fn HardApproachLaw::fit()` on a strictly decreasing linear neg-log grid returns a power-law exponent `b > 0` instead of the linear case `b = 0`, so the linear refit that would recover the `slope` never runs and the boundary shape is misrepresented.

## Symptom and reproduction

No current test exercises a falling (negative-slope) line: the `fit` tests in [`packages/treetime-grid/src/__tests__/hard_approach_law.rs`](../../packages/treetime-grid/src/__tests__/hard_approach_law.rs) cover rising power laws (`test_hard_approach_law_fit_recovers_exponent`) and rising lines only (`test_hard_approach_law_fit_finite_line_recovers_slope` with `slope = 0.5`, `test_hard_approach_law_fit_clamps_negative_exponent_to_zero` with `slope = 2.0`). The defect is therefore latent, not a live failure.

To reproduce, fit a grid storing the exact line `y(t) = -2 t + 4` over `t ∈ [0.1, 1.0]` with `t_hard = 0`. The correct recovery is `b = 0`, `slope = -2.0`; the fit instead returns a spurious `b > 0` with `slope = 0`, because the stage-1 log-distance regression slope is negative for a falling line (see root cause). Adding this case to the `fit` tests is the acceptance test for the fix.

## Root cause

`fit` classifies the boundary shape in stage 1 by least-squares regression of `(ln|t - t_hard|, y)`, taking the slope as `-b` [`packages/treetime-grid/src/hard_approach_law.rs#L84-L101`](../../packages/treetime-grid/src/hard_approach_law.rs#L84-L101). Over a fit window bounded away from `t_hard`, a falling line `y` decreases as `ln|t - t_hard|` increases, so the regression slope is negative and `b` comes out spuriously positive (here `≈ 0.516`). Because the linear refit only runs when `b <= B_THRESHOLD` (`0.01`) [`packages/treetime-grid/src/hard_approach_law.rs#L103-L128`](../../packages/treetime-grid/src/hard_approach_law.rs#L103-L128), a genuinely linear but falling density never reaches it and is reported as a power law. A rising line gives a negative raw `b` that clamps to `0` and passes; a falling line does not.

## Impact

Latent today. A production path now fits a `HardApproachLaw`: `create_branch_length_likelihood` fits one on the left boundary at `t_hard = 0` [`packages/treetime/src/timetree/inference/branch_length_likelihood.rs`](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs). For a zero-mutation branch (`n = 0`) the neg-log density there is `y = mu*t + C`, a rising line, so the falling case does not arise on that path and no current output is wrong. The defect matters when a falling neg-log line reaches the fit (for example a composed or convolved message whose finite side decreases toward the boundary): misclassifying it as a vanishing power law would misplace the boundary approach and any mode sitting on it.

## Fix approach

The `b >= 0` clamp cannot distinguish a falling line (spurious positive `b`) from a true power law. Detect linearity independently of the log-distance fit, for example by comparing the residuals of a linear-in-`t` fit against the linear-in-`ln|dt|` fit and selecting the better model, or by restricting the power-law branch to windows where `y` increases toward `t_hard`.
