# Skyline dates and Tc drift over refinement rounds at low signal

With `--coalescent skyline` the node dates and the fitted $T_c$ degrade with every refinement round at low signal, also when the clock rate is fixed at its true value. `--coalescent opt` does not show it. The default of two rounds (`--max-iter 2`) hides most of it.

## Evidence

Clean validation data (`validation/data`, 400 tips, 10 replicates of the `const`, `exp` and `piecewise` schedules), $L=10^3$, true topology, `--keep-root --clock-filter 0`; median over runs:

| variant | $\mu$ | rounds | date RMSE (years) | 90% quantile | $\ln(T_c/T_c^\text{true})$ |
| --- | --- | --- | --- | --- | --- |
| `keep-skyline-fixed` (true rate), 30 runs | $3\cdot10^{-4}$ | 2 | 0.24 | 0.52 | $-0.11$ |
| | | 6 | 2.77 | 6.02 | $+1.93$ |
| `keep-skyline-est`, 30 runs | $3\cdot10^{-4}$ | 2 / 6 | 0.24 / 1.92 | | |
| `keep-skyline-est` | $10^{-3}$ | 2 / 6 | 0.20 / 0.39 | | |
| `keep-opt-est` | $3\cdot10^{-4}$ | 2 / 6 | 0.33 / 0.34 | | |

The estimated-rate rows are for the shrinking rate window ([kb/decisions/timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md#shrinking-probe-window)); with the fixed window they are 0.24 / 2.04 and 0.20 / 0.54. At $\mu \ge 3\cdot10^{-3}$ six rounds change the date RMSE by less than 0.04 years.

## Locations

- The skyline $T_c(t)$ is refitted every round from the expected lineage counts of the current posteriors (`fn round_models()` calls `fn coalescent_timescale()`, [packages/treetime/src/timetree/refinement_loop.rs#L118](../../packages/treetime/src/timetree/refinement_loop.rs#L118); `fn optimize_skyline()` [packages/treetime/src/coalescent/skyline.rs#L15](../../packages/treetime/src/coalescent/skyline.rs#L15)), while the prior holds the expected lineage counts frozen after the first rebuild ([kb/decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md))
- The skyline segments are equal-width over the span of the expected lineage counts ([packages/treetime/src/coalescent/skyline.rs#L40-L41](../../packages/treetime/src/coalescent/skyline.rs#L40-L41)), so the segments move with the dates

> [!IMPORTANT]
> **Investigation required.** The mechanism. $T_c$ grows about seven-fold over six rounds at the true rate while the dates spread, which points to a feedback between the refitted skyline and the node times that the single $T_c$ of `opt` does not have. Trace $T_c(t)$, the segment boundaries and the root date per round on a few `const` replicates at $L=10^3$, $\mu=3\cdot10^{-4}$ with `--clock-rate 3e-4`, and repeat with the segment boundaries fixed after the first round.

## Validation

`keep-skyline-fixed` and `keep-skyline-est` at $L=10^3$, $\mu \in \{3\cdot10^{-4}, 10^{-3}\}$ with `--max-iter 6`: the date RMSE and $T_c$ bias after six rounds no larger than after two.
