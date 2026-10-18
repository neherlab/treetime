# Coalescent statistics come from two different lineage counts

Phase 2 of [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md).

## Problem

The prior and the $T_c$ fit use the expected lineage counts (`LineageExpectation`); the polytomy merger rates, the
convergence metric `log_lh_coal`, the reported coalescent likelihood, the skyline segment layout and the $T_c$ band use
point-time lineage counts in `CoalescentModel`, with the root patch `LINEAGES_BEFORE_ROOT` and the per-lineage rate
floor `max(0.5, k - 1)`.

## Fix direction

Take all of them from the expectation: merger rates $(\bar k - 1)/(2T_c)$, the expected complete-data
log-likelihood $-\sum_i I_i/T_{c,i} - \sum_i M_i \ln T_{c,i}$ for reporting and convergence, the band from the Hessian
of the expected objective, the segment range from the expectation grid. Remove `CoalescentModel`, the event and
lineage-count code and the edge-based total likelihood.

## Validation

Tests against analytic expected statistics; outputs documented as a v0 divergence.
