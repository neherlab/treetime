# v0/v1 Deviations

Deliberate behavioral differences between v1 (Rust) and v0 (Python), with rationale. One file per deviation.

Distinct from [known issues](../port-known-issues/_index.md) (unintentional differences, not-yet-ported features) and [proposals](../port-proposals/_index.md) (new v1 features not in v0, pre-implementation).

## Summary

| Domain     | Deviation                                                                                                                                                       | Impact                                               |
| ---------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------- |
| Ancestral  | [Joint ML reconstruction removed](ancestral-joint-reconstruction-removed.md#joint-ml-ancestral-reconstruction-removed)                                          | `--method-anc joint` panics                          |
| Ancestral  | [Fitch root ambiguity: deterministic selection](ancestral-fitch-deterministic-root-state.md#fitch-root-ambiguity-deterministic-selection)                       | Root sequence at tied positions                      |
| GTR        | [Uninformative root_state filtering](gtr-uninformative-root-state-filtering.md#uninformative-root_state-filtering-in-gtr-inference)                             | pi, W, mu shift for gappy datasets                   |
| Coalescent | [Coalescent multiplication ordering](coalescent-multiplication-ordering.md#coalescent-multiplication-ordering)                                                  | Grid points, interpolation path                      |
| Sequences  | [Dense and sparse sequence representation](sequence-representation-dense-sparse.md#dense-and-sparse-sequence-representation)                                    | Memory usage, large alignments                       |
| Structure  | [Graph-based phylogenetic representation](graph-based-phylogenetic-representation.md#graph-based-phylogenetic-representation)                                   | Tree storage, traversal, data access                 |
| Optimize   | [6-method branch optimization selection](optimize-newton-raphson-per-edge.md#per-edge-branch-length-optimization-6-method-selection)                            | Per-branch convergence speed, method choice          |
| Structure  | [Partition system architecture](partition-system-architecture.md#partition-system-architecture)                                                                 | Code organization, extensibility                     |
| Commands   | [Standalone branch length optimization command](command-optimize-standalone.md#standalone-branch-length-optimization-command)                                   | New `optimize` subcommand                            |
| Commands   | [Standalone tree pruning command](command-prune-standalone.md#standalone-tree-pruning-command)                                                                  | New `prune` subcommand                               |
| I/O        | [Multi-format tree I/O](multi-format-tree-io.md#multi-format-tree-io)                                                                                           | 8 formats, interop with UShER/Auspice/PhyloXML       |
| Clock      | [Pre-filter clock allows negative rate](clock-prefilter-relaxed-positive-rate.md#clock-pre-filter-allows-negative-rate-during-root-finding)                     | Outlier detection proceeds on negative-rate datasets |
| Coalescent | [Actual multiplicity instead of fixed 2.0](coalescent-total-lh-actual-multiplicity.md#coalescent-total-lh-uses-actual-multiplicity-instead-of-fixed-2)          | Correct merger rate for polytomies                   |
| GTR        | [Sparse fixed-position scalar rate approximation](sparse-fixed-position-scalar-rate-approximation.md#sparse-fixed-position-scalar-rate-approximation)           | Fixed positions use scalar mu, not per-site rate     |
| Mugration  | [Pseudo-count smoothing on initial pi](mugration-pseudo-count-initial-pi.md#pseudo-count-smoothing-on-initial-mugration-pi)                                     | Smoother prior for first reconstruction              |
| Mugration  | [Root state uniform-threshold filtering](mugration-root-state-filtering.md#root-state-uniform-threshold-filtering-in-mugration-gtr-inference)                   | Avoids state-order bias at uninformative root        |
| Optimize   | [Dense initial guess uses discrete subs](optimize-dense-initial-guess-soft-hamming.md#dense-initial-branch-length-guess-uses-soft-hamming-on-per-edge-messages) | Initial branch length from argmax substitutions      |
| Optimize   | [Indel contribution to branch length likelihood](optimize-indel-contribution-to-likelihood.md#indel-contribution-to-branch-length-likelihood)                   | Positive branch length for indel-only edges          |
| Timetree   | [CI output not gated on --confidence](timetree-ci-output-ungated.md#confidence-interval-output-not-gated-on---confidence)                                       | Extra CI file when using marginal modes              |
| Optimize   | [Model-aware zero-branch shortcut](optimize-model-aware-zero-branch-shortcut.md#model-aware-zero-branch-shortcut)                                               | Non-unimodal models bypass derivative shortcut       |
| Timetree   | [Timetree EM loop does not collapse zero-length branches](timetree-no-zero-branch-collapse-in-loop.md#timetree-em-loop-does-not-collapse-zero-length-branches)  | Topology cleanup in optimize pre-step, not EM loop   |
| Prune      | [Merged-sibling branch length uses Jukes-Cantor correction](prune-merge-jukes-cantor-branch-length.md#merged-sibling-branch-length-uses-jukes-cantor-correction-not-raw-p-distance) | Corrected evolutionary distance on shared-mutation edges |
