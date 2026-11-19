# The Tc confidence band under-covers on simulated data

> [!IMPORTANT]
> **Decision required.** The approved band [kb/decisions/coalescent-skyline-hessian-confidence-bands.md](../decisions/coalescent-skyline-hessian-confidence-bands.md) is a local Gaussian approximation with node times held fixed, and it states no coverage guarantee. The options:
>
> - Keep the band and document the limit in the `--help` text and the output description. No extra cost, but the band stays somewhat too narrow on short alignments (0.90 on the `const` schedule)
> - Include node-time uncertainty, for example by refitting $T_c$ over node times sampled from their marginal posteriors. Coverage improves, at the cost of extra passes per run
>
> Evidence: the coverage below.

## Problem

With `--coalescent opt` and the default width of 2 standard deviations of $\ln T_c$ [packages/app-commands/src/commands/timetree/args.rs#L885-L896](../../packages/app-commands/src/commands/timetree/args.rs#L885-L896), about 95% under the Gaussian approximation, the band covers the simulated coalescent time scale $T_c$ in 63% of the runs at alignment length $L = 1000$. The numbers come from the timetree suite of [validation/README.md](../../validation/README.md), with the root kept and the true clock rate (column `tc_coverage`).

The 63% pools the `exp` and `piecewise` schedules, whose true $T_c(t)$ varies in time, so a constant-$T_c$ band cannot cover it there. Rerun on `rust` of 2026-10-24 (`keep-opt-fixed`, $L = 1000$, $\mu \in \{10^{-3}, 3 \cdot 10^{-3}\}$, rep00 to rep09, 60 runs): on the `const` schedule, where the model is right, the band covers in 18 of 20 runs (0.90 against about 0.95); on `exp` in none (log bias +0.56) and on `piecewise` in 0.6 to 0.9 (log bias +0.2 to +0.84). $\mu = 3 \cdot 10^{-4}$ is not rerun. The remaining question is the modest undercoverage on well-specified data and what the `--help` text states.

The band is computed from the curvature of the coalescent likelihood in $\ln T_c$ with the node times held fixed (`fn skyline_confidence_band()` [packages/treetime/src/coalescent/skyline.rs#L250-L261](../../packages/treetime/src/coalescent/skyline.rs#L250-L261)). The `--help` text of `--coalescent-confidence` describes this curvature band, but it does not state that the uncertainty of the node times is left out.

## Validation

Depends on the decision: the documented limit, or coverage close to the nominal level on the validation suite.

## Proposed implementation plan

The locations are current (`coalescent_confidence` with its `--help` text at `args.rs` L885-L896, `fn skyline_confidence_band()` at `skyline.rs` L254). The help text still does not state that node-time uncertainty is left out. The plan depends on the open decision.

- **Keep the band**: extend the doc comment of `coalescent_confidence` in `packages/app-commands/src/commands/timetree/args.rs` (and of the `clock` command if it shares the flag) to state that the band holds the node times fixed and so under-covers; run `just gen` and `just generated-check`; no test, because the change is help text only
- **Include node-time uncertainty**: record the approved method in `kb/decisions/coalescent-skyline-hessian-confidence-bands.md`; sample node times from their marginal posteriors with a seeded `StdRng`, refit $\ln T_c$ per sample, and combine the curvature variance with the variance across samples; tests on a simulated `const` tree from `treetime-simulation` in the slow tier (`test_sim_`) asserting coverage near 0.95 over many seeded replicates, with the sample size and false-fail rate stated (oracle: the simulated $T_c$)
- Validation for both: `just check`, `just test-rs`; for the second option, the `keep-opt-fixed` rerun of the validation suite at $L = 1000$
