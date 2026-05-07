# Proposals

New v1 features not present in v0. Each proposal records scientific motivation, expected impact, related code, and a validation plan before implementation. Once implemented, the entry moves to [`decisions/`](../decisions/README.md).

v0 features not yet implemented belong in [`issues/`](../issues/README.md) (all types) and [`algo/unimplemented.md`](../algo/unimplemented.md) (algorithms only), not here.

## Summary

| Domain    | Proposal                                                                                                                                                         | Motivation                                                        |
| --------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------- |
| Ancestral | [Iterative outer GTR refinement for ancestral reconstruction](ancestral-iterative-gtr-refinement.md#iterative-outer-gtr-refinement-for-ancestral-reconstruction) | Improve latent-state and model co-estimation when `--model infer` |
| Mugration | [Full forward-backward reconstruction per iteration](mugration-full-reconstruction-per-iteration.md)                                                             | Correct EM: refresh both message directions each iteration        |
| Optimize  | [Convergence and robustness](optimize-convergence-and-robustness.md)                                                                                             | Sparse fix, damping, acceleration, per-edge safety, GTR inference |
| Timetree  | [Propagate convergence metric evaluation errors](convergence-metric-error-propagation.md)                                                                        | Prevent silent exclusion of failed likelihood components          |
| GTR       | [Configurable convergence norm for GTR inference](gtr-inference-convergence-norm.md)                                                                             | Dimension-independent convergence for site-specific inference     |
| GTR       | [Max-rate interpolation grid for site-specific GTR](gtr-site-specific-interpolation-grid.md)                                                                     | Better interpolation accuracy for heterogeneous rate sites        |
| I/O       | [Configuration file format for multi-partition analysis](config-file-multi-partition.md)                                                                         | Multi-segment genomes, mixed data types, per-partition models     |
| I/O       | [Unified input format support for analysis commands](unified-input-format-support.md)                                                                            | Bypass name matching via Auspice/MAT/PhyloXML direct input        |
| Optimize  | ~~Model-aware zero-branch shortcut~~ [implemented](../decisions/optimize-model-aware-zero-branch-shortcut.md)                                                    | Restrict derivative shortcut to unimodal models                   |
| Optimize  | [Indel model alternatives](optimize-indel-model-alternatives.md)                                                                                                 | Length-weighted, separate ins/del, affine-inspired indel models   |
| Optimize  | [Per-iteration indel rate re-estimation](optimize-indel-rate-per-iteration-reestimation.md)                                                                      | Restore ECM coordinate-ascent re-estimation of indel rate         |
| Graph     | [Pure transform architecture for graph traversals](graph-pure-transform-architecture.md#pure-transform-architecture-for-graph-traversals)                        | Explicit data flow, no locks, type-safe pipeline staging          |
