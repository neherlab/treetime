# Timetree loop starts from a prior-free pass that is too old at low signal

Phase 1 of [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md).

## Problem

`fn run_initial_round()` ([`round.rs`](../../packages/treetime/src/timetree/round.rs)) runs the time passes without a
prior, fits $T_c$ and the expected lineage counts to them, runs a prior pass and rebuilds. The prior-free node times
are too old at low signal even at the true rate (validation, $L=10^3$, $\mu=3\cdot10^{-4}$: bias $-2.7$ years, RMSE
3.65 against 0.47 for Gaussian marginals at the regression rate; lineage-curve distance 1043 against 65
lineage-years).

## Fix direction

Merge `claude-hard-cell-interp`. Seed the expectation from the Gaussian node marginals of 0b (truncated at each node's
earliest tip, the root possibly at the regression's root time) and the first $T_c$ by the EM fit on it; the initial
round becomes round 1 of the loop, with the expectation rebuilt once from its posteriors and then frozen.

## Validation

Calibration and accuracy equal or better than mean field + EM; one time pass fewer; `data/dengue/500` completes.
Closes [M-timetree-initial-round-applies-coalescent-prior.md](M-timetree-initial-round-applies-coalescent-prior.md).
