# Port Proposals

New v1 features that v0 does not have. Each proposal records scientific motivation, expected impact, related code, and a validation plan before implementation. Once implemented, the deviation from v0 moves to [`docs/port-intentional-changes/`](../port-intentional-changes/_index.md).

v0 features not yet ported belong in [`docs/port-known-issues/`](../port-known-issues/_index.md) (all types) and [`docs/port-algo-inventory/unimplemented.md`](../port-algo-inventory/unimplemented.md) (algorithms only), not here.

## Summary

| Domain    | Proposal                                                                                                                                                                       | Motivation                                                        |
| --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------- |
| Ancestral | [Iterative outer GTR refinement for ancestral reconstruction](ancestral-iterative-gtr-refinement.md#iterative-outer-gtr-refinement-for-ancestral-reconstruction)               | Improve latent-state and model co-estimation when `--model infer` |
| Mugration | [Full forward-backward reconstruction per iteration](mugration-full-reconstruction-per-iteration.md)                                                                           | Correct EM: refresh both message directions each iteration        |
| Optimize  | [Reliable convergence and method choice](optimize-convergence-and-method-choice.md#reliable-convergence-and-method-choice-for-branch-length-optimization)                      | Damping, GTR wiring, gap-aware init, Brent alternative            |
| Timetree  | [Propagate convergence metric evaluation errors](convergence-metric-error-propagation.md)                                                                                      | Prevent silent exclusion of failed likelihood components          |
| GTR       | [Configurable convergence norm for GTR inference](gtr-inference-convergence-norm.md)                                                                                           | Dimension-independent convergence for site-specific inference     |
| GTR       | [Max-rate interpolation grid for site-specific GTR](gtr-site-specific-interpolation-grid.md)                                                                                   | Better interpolation accuracy for heterogeneous rate sites        |
| CLI       | [Configuration file format for multi-partition analysis](config-file-multi-partition.md)                                                                                       | Multi-segment genomes, mixed data types, per-partition models     |
| Optimize  | ~~[Model-aware zero-branch shortcut](optimize-model-aware-zero-branch-shortcut.md)~~ → [implemented](../port-intentional-changes/optimize-model-aware-zero-branch-shortcut.md) | Restrict derivative shortcut to unimodal models                   |
