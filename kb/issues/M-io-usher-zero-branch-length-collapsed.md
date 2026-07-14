# UShER input collapses explicit zero branch length into absence

The UShER MAT reader converts an explicit branch length of zero to `None`. Zero substitutions and an unspecified branch length are distinct states; collapsing them prevents exact round trips and can change downstream defaulting behavior.

The embedded Newick parser produces an optional branch length, but `fn usher_to_graph()` converts absence to `0.0` in `UsherNodeImpl` [packages/treetime-io/src/usher_mat.rs#L232-L247](../../packages/treetime-io/src/usher_mat.rs#L232-L247). `fn TreeIrUsherReader::usher_node_to_graph_components()` then maps every zero back to `None` [packages/treetime-io/src/tree_ir/usher.rs#L117-L142](../../packages/treetime-io/src/tree_ir/usher.rs#L117-L142).

## Potential solutions

- O1. Preserve presence independently of numeric value.
- O2. Model absent/zero/positive branch lengths as an explicit enum. This is useful only if the distinction must propagate beyond existing `Option<f64>` boundaries.

## Recommendation

Preserve protobuf presence independently of numeric value. `Some(0.0)` remains `Some(0.0)`, while an absent field remains `None`. Reject non-finite or negative values according to the common branch-length contract.

## Related issues

- [M-io-usher-mat-mutation-loss-is-implicit.md](M-io-usher-mat-mutation-loss-is-implicit.md)
- [N-core-branch-length-clamping.md](N-core-branch-length-clamping.md)
