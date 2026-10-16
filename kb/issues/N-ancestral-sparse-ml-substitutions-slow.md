# Sparse ML substitution extraction sorts positions per edge

`fn compute_ml_subs_for_nodes()` [packages/treetime/src/partition/marginal/sparse/forward.rs#L228-L254](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L228-L254) finds the substitutions between the reconstructed parent and child sequences of one edge. It chains the positions of the Fitch substitutions of the edge with the explicit positions of both node profiles, sorts and deduplicates them, and calls `fn map_state()` twice per position, once for the parent and once for the child. Each call searches the explicit positions of a node profile.

On `lassa/L/500`, sparse `ancestral` with one thread spends 10.8 % of all CPU samples in this function.

> [!IMPORTANT]
> **Investigation required.** The explicit positions of both profiles are sorted by position, so a merge of the sorted lists could replace the sort, and the same walk could read both node states without the searches of `map_state()`. The gain is not measured.
