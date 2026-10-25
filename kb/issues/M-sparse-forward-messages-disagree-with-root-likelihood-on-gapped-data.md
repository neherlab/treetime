# Sparse forward messages disagree with the root likelihood on gapped data

> [!IMPORTANT]
> **Investigation required.** The deviations were measured with debug code in `build/likelihood.rs` that is not committed. Before a fix, commit a test that checks the identity below at every node of a small tree with gaps and unknown leaves, against the dense backend as the oracle.

## Summary

For a correct marginal reconstruction, combining the outside message at a node with the subtree messages of its children gives the full-tree likelihood at every node. In the sparse backend this holds on data without gaps and fails on data with gaps, deletions, and unknown stretches. The node log-likelihoods differ from the root log-likelihood by up to:

| Dataset | Max \|node log_lh − tree log_lh\| |
| --- | --- |
| `flu/h3n2/500` | 0.002 |
| `rsv/a/2000` | 50 |
| `dengue/2000` | 583 |

The node value was computed as `combine_messages()` of `forward[e].msg_to_child` for the parent edge and `backward[c].msg_from_child` for each child edge. The tree value is the root log-likelihood of the same `PartitionMarginalSparse` after `marginal_update()`, with GTR and `IndelPolicy::Exclude`. The deviation grows with the number of gapped and unknown positions in the dataset.

## Effect

Any consumer of the forward messages takes the wrong likelihood: branch-length optimization uses the messages at the two ends of each edge, and so would a likelihood placement or SPR search. The likelihood SPR prototype (`build/likelihood.rs`) works around this with a constant shift per edge and per node, so that each local value matches the root value at the current topology. With the shifts, the reattachment error is 6 to 8 log units on rsv and dengue and 27 to 32 on Lassa; on flu it is 1e-12. That error is larger than the average predicted gain of a move, so the search cannot rank moves reliably on gapped data.

## Suspected cause

Two mechanisms are likely:

- The backward pass masks the non-character positions of a child (gap or unknown), through `fill_missing_references()`. The forward pass may treat those positions differently, so the two directions disagree on which sites carry information on an edge.
- Resolved sites are dropped from the variable set below `EPS = 1e-4`. This makes the combine step non-associative, but by an amount that does not explain deviations of hundreds of log units.

[Sparse marginal likelihood depends on the branch length of an unknown leaf](H-sparse-marginal-likelihood-depends-on-branch-of-unknown-leaf.md) is probably the same defect in a four-taxon case.

## Related

- [Sparse root invariance violation](M-ancestral-sparse-root-invariance.md)
- [Sparse marginal cavity compression unverified](N-sparse-marginal-cavity-compression-unverified.md)
