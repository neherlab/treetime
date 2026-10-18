# Timetree loop: a non-positive regression rate stops the run before the init phase

Phase 1 of [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md).
The init phase itself is done (`fn initialize()`, [`init.rs`](../../packages/treetime/src/timetree/init.rs)): it seeds
the expected lineage counts from the Gaussian node marginals at the regression rate, fits $T_c$ to them, infers the
node times once under that prior and rebuilds the expectation, which the loop keeps frozen; results in the proposal,
"Phase 1 result".

## Problem

When the root-to-tip regression gives a non-positive rate and no `--clock-rate` is given, timetree stops before the
init phase. Four runs of the phase 0b and phase 1 gates (`const` rep03 and rep09, $L=10^3$, `--keep-root`) stop this
way, with every variant of the init.

## Fix direction

Start from the Gaussian rate estimate ($\hat\mu$ from the $D{:}M{:}N$ split) instead of stopping: in the evaluation of
the proposal it never failed and was less noisy than OLS, though 0.4-26% high, and the marginal-likelihood step
corrects the rate from there. The regression stays the start whenever its rate is positive.

## Validation

The four runs complete with a rate error comparable to their neighbors; no other run changes.
