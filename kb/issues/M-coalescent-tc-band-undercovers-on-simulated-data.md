# The Tc confidence band under-covers on simulated data

> [!IMPORTANT]
> **Decision required.** The approved band [kb/decisions/coalescent-skyline-hessian-confidence-bands.md](../decisions/coalescent-skyline-hessian-confidence-bands.md) is a local Gaussian approximation with node times held fixed, and it states no coverage guarantee. The options:
>
> - Keep the band and document the limit in the `--help` text and the output description. No extra cost, but the band stays too narrow for short alignments
> - Include node-time uncertainty, for example by refitting $T_c$ over node times sampled from their marginal posteriors. Coverage improves, at the cost of extra passes per run
>
> Evidence: the coverage below.

## Problem

With `--coalescent-opt` and the default width of 2 standard deviations of $\ln T_c$ [packages/app-commands/src/commands/timetree/args.rs#L426-L428](../../packages/app-commands/src/commands/timetree/args.rs#L426-L428), about 95% under the Gaussian approximation, the band covers the simulated coalescent time scale $T_c$ in 63% of the runs at alignment length $L = 1000$. The numbers come from the timetree suite of [validation/README.md](../../validation/README.md), with the root kept and the true clock rate (column `tc_coverage`).

## Validation

Depends on the decision: the documented limit, or coverage close to the nominal level on the validation suite.

## Related issues

- [M-coalescent-skyline-underestimates-tc-on-simulated-data.md](M-coalescent-skyline-underestimates-tc-on-simulated-data.md)
