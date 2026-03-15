# Port Proposals

Candidate behavioral changes or algorithmic extensions for v1 that are grounded in the current codebase but are not yet accepted or implemented.

These documents are the proposal-stage analogue of [`docs/port-intentional-changes/`](../port-intentional-changes/_index.md). They record scientific motivation, current behavior, related code, expected impact, and validation ideas before any decision is made.

## Summary

| Domain    | Proposal                                                                                                                                                         | Motivation                                                        |
| --------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------- |
| Ancestral | [Iterative outer GTR refinement for ancestral reconstruction](ancestral-iterative-gtr-refinement.md#iterative-outer-gtr-refinement-for-ancestral-reconstruction) | Improve latent-state and model co-estimation when `--model infer` |
| Optimize  | [Reliable convergence and method choice](optimize-convergence-and-method-choice.md#reliable-convergence-and-method-choice-for-branch-length-optimization)        | Damping, GTR wiring, gap-aware init, Brent alternative            |
