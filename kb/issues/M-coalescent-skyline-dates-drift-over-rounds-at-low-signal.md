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

- The skyline $T_c(t)$ is refitted every round from the expected lineage counts of the current posteriors (`fn round_models()` calls `fn coalescent_timescale()`, [packages/treetime/src/timetree/refinement_loop.rs#L107](../../packages/treetime/src/timetree/refinement_loop.rs#L107); `fn optimize_skyline()` [packages/treetime/src/coalescent/skyline.rs#L15](../../packages/treetime/src/coalescent/skyline.rs#L15)), while the prior holds the expected lineage counts frozen after the first rebuild ([kb/decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md))
- The skyline segments are equal-width over the span of the expected lineage counts ([packages/treetime/src/coalescent/skyline.rs#L40-L41](../../packages/treetime/src/coalescent/skyline.rs#L40-L41)), so the segments move with the dates

> [!IMPORTANT]
> **Investigation required.** The mechanism. $T_c$ grows about seven-fold over six rounds at the true rate while the dates spread, which points to a feedback between the refitted skyline and the node times that the single $T_c$ of `opt` does not have. Trace $T_c(t)$, the segment boundaries and the root date per round on a few `const` replicates at $L=10^3$, $\mu=3\cdot10^{-4}$ with `--clock-rate 3e-4`, and repeat with the segment boundaries fixed after the first round.

## Validation

`keep-skyline-fixed` and `keep-skyline-est` at $L=10^3$, $\mu \in \{3\cdot10^{-4}, 10^{-3}\}$ with `--max-iter 6`: the date RMSE and $T_c$ bias after six rounds no larger than after two.

## Proposed implementation plan

The locations are current (`refinement_loop.rs` L107 refits the timescale each round when `mode.is_optimized()`; `fn optimize_skyline()` at `skyline.rs` L15, equal-width boundaries at L40-L41). The fix depends on the mechanism of the Investigation required block, so the plan starts with it.

1. Trace per round, on 3 `const` replicates at $L=10^3$, $\mu=3\cdot10^{-4}$ with `--clock-rate 3e-4 --max-iter 6`: the segment boundaries, $T_c(t)$ per segment, the expected lineage counts that the fit reads, and the root date (temporary `debug!` output in `fn round_models()`, removed after the trace)
2. Repeat with the segment boundaries fixed after the first round (temporary change in `fn optimize_skyline()`), and with the lineage counts of the fit frozen as the prior's are, to tell which feedback drives the drift
3. Fix the mechanism that the trace shows. Candidates: fix the segment boundaries after the first fit; fit the skyline on the same frozen lineage counts as the prior (`kb/decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md`). A change of the skyline estimator is a scientific choice and needs user approval and a `kb/decisions/` entry
4. Tests: a regression test in `packages/treetime/src/coalescent/__tests__/` on a small simulated tree from `treetime-simulation` (constant $T_c$, low signal, fixed seed) asserting that $T_c$ after six rounds lies no further from the true value than after two (oracle: the known simulated $T_c$); add it to the slow tier (`test_sim_`) when it needs replicates
5. Validation: the Validation section above, plus `just check` and `just test-rs`
