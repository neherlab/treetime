# Port Proposals

New v1 features that v0 does not have. Each proposal records scientific motivation, expected impact, related code, and a validation plan before implementation. Once implemented, the deviation from v0 moves to [`docs/port-intentional-changes/`](../port-intentional-changes/_index.md).

v0 features not yet ported belong in [`docs/port-known-issues/`](../port-known-issues/_index.md) (all types) and [`docs/port-algo-inventory/unimplemented.md`](../port-algo-inventory/unimplemented.md) (algorithms only), not here.

## Summary

| Domain    | Proposal                                                                                                                                                         | Motivation                                                        |
| --------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------- |
| Ancestral | [Iterative outer GTR refinement for ancestral reconstruction](ancestral-iterative-gtr-refinement.md#iterative-outer-gtr-refinement-for-ancestral-reconstruction) | Improve latent-state and model co-estimation when `--model infer` |
| Optimize  | [Reliable convergence and method choice](optimize-convergence-and-method-choice.md#reliable-convergence-and-method-choice-for-branch-length-optimization)        | Damping, GTR wiring, gap-aware init, Brent alternative            |
| Timetree  | [Convergence metric error propagation](convergence-metric-error-propagation.md#propagate-convergence-metric-evaluation-errors-instead-of-swallowing-them)        | Failed metric components silently excluded from lh_total          |
