# v0/v1 Deviations

Deliberate behavioral differences between v1 (Rust) and v0 (Python), with rationale. One file per deviation.

Distinct from [known issues](../port-known-issues/_index.md) (unintentional differences, not-yet-ported features) and [proposals](../port-proposals/_index.md) (new v1 features not in v0, pre-implementation).

## Summary

| Domain     | Deviation                                                                                                                                                | Impact                               |
| ---------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------ |
| Ancestral  | [Joint ML reconstruction removed](ancestral-joint-reconstruction-removed.md#joint-ml-ancestral-reconstruction-removed)                                   | `--method-anc joint` panics          |
| Ancestral  | [Fitch root ambiguity: deterministic selection](ancestral-fitch-deterministic-root-state.md#fitch-root-ambiguity-deterministic-selection)                | Root sequence at tied positions      |
| GTR        | [Uninformative root_state filtering](gtr-uninformative-root-state-filtering.md#uninformative-root_state-filtering)                                       | pi, W, mu shift for gappy datasets   |
| Coalescent | [Coalescent multiplication ordering](coalescent-multiplication-ordering.md#coalescent-multiplication-ordering)                                           | Grid points, interpolation path      |
| Sequences  | [Dense and sparse sequence representation](sequence-representation-dense-sparse.md#dense-and-sparse-sequence-representation)                             | Memory usage, large alignments       |
| Structure  | [Graph-based tree structure](tree-structure-graph-based.md#graph-based-tree-structure)                                                                   | Tree storage, traversal, data access |
| Optimize   | [Newton-Raphson per-edge instead of Brent](optimize-newton-raphson-per-edge.md#per-edge-branch-length-optimization-uses-newton-raphson-instead-of-brent) | Per-branch convergence speed         |
| Structure  | [Partition system architecture](partition-system-architecture.md#partition-system-architecture)                                                          | Code organization, extensibility     |
