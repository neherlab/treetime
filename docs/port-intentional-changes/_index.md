# v0/v1 Deviations

Intentional behavioral differences between v1 (Rust) and v0 (Python).

Deviations are distinct from known issues (unintentional differences to be addressed) and unimplemented features (capability not yet ported).

## Summary

| Domain     | Deviation                                                                                                                                 | Impact                               |
| ---------- | ----------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------ |
| Ancestral  | [Joint ML reconstruction removed](ancestral-joint-reconstruction-removed.md#joint-ml-ancestral-reconstruction-removed)                    | `--method-anc joint` panics          |
| Ancestral  | [Fitch root ambiguity: deterministic selection](ancestral-fitch-deterministic-root-state.md#fitch-root-ambiguity-deterministic-selection) | Root sequence at tied positions      |
| GTR        | [Uninformative root_state filtering](gtr-uninformative-root-state-filtering.md#uninformative-root_state-filtering)                        | pi, W, mu shift for gappy datasets   |
| Coalescent | [Coalescent multiplication ordering](coalescent-multiplication-ordering.md#coalescent-multiplication-ordering)                            | Grid points, interpolation path      |
| Sequences  | [Dense and sparse sequence representation](sequence-representation-dense-sparse.md#dense-and-sparse-sequence-representation)              | Memory usage, large alignments       |
| Structure  | [Graph-based tree structure](tree-structure-graph-based.md#graph-based-tree-structure)                                                    | Tree storage, traversal, data access |
| Structure  | [Partition system architecture](partition-system-architecture.md#partition-system-architecture)                                           | Code organization, extensibility     |
