# Sparse posterior of a non-root node folds sites without data below a child

At a position where one child subtree of an internal node has no data, the sparse posterior of the node differs from the dense posterior when the node is not the root. The forward pass folds such a position into the shared profile of its state once the posterior exceeds $1 - 10^{-4}$, although the position carries a different posterior than the shared profile. The root keeps the same positions explicit and agrees with dense.

## Mechanism

The backward pass keeps a child-missing site of a node explicit in the node's message to its parent: the missing child has a gap or unknown reference in the candidate site table of `fn combine_child_messages()` ([packages/treetime/src/partition/marginal/sparse/backward.rs](../../packages/treetime/src/partition/marginal/sparse/backward.rs)), so `fn combine_messages()` ([packages/treetime/src/partition/marginal/sparse/message.rs](../../packages/treetime/src/partition/marginal/sparse/message.rs)) sees a reference that differs from the node state. At the root, this product is the posterior, so the position stays explicit.

Below the root, `fn SparseForwardPass::edge()` ([packages/treetime/src/partition/marginal/sparse/forward.rs](../../packages/treetime/src/partition/marginal/sparse/forward.rs)) rebuilds the posterior from two messages, the message from the parent and the node's own message to the parent. The candidate site table of this product carries the label of the node's explicit entry and no reference of the missing child, so every message agrees with the node state. A position whose posterior exceeds $1 - 10^{-4}$ then folds into the shared profile of its state, which assumes that every child has data there. `fn compute_msg_to_child()` divides this shared profile by the message of each child, so the outside messages of all children of the node, and their branch-length objectives, lose the difference.

[kb/algo/ancestral.md](../algo/ancestral.md) states that a child-missing site stays explicit at the node and folds only at an ancestor. The forward pass does not follow this rule for the posterior of a non-root node.

## Evidence

The dense representation is the oracle. `test_root_polytomy_sparse_matches_dense` in [packages/treetime/src/partition/marginal/sparse/__tests__/test_root_polytomy.rs](../../packages/treetime/src/partition/marginal/sparse/__tests__/test_root_polytomy.rs) compares the root log-likelihood, the internal-node posteriors, the change of every edge objective at three branch lengths, and the optimal branch length of every edge with dense, for a polytomy of 8 or 12 tips with runs of `N` at the root and one level below it. The case `long_branches_below_root` (8 tips at branch length 0.8, one run of `N` per tip except one, polytomy node `P` at 0.02 below the root) is ignored because of this issue:

| Quantity | Largest difference, sparse minus dense |
| --- | ---: |
| Edge objective change at t = 0.001 (edge `c6`) | -0.043 nats |
| Optimal branch length (edge `c6`) | 4.4e-5 |
| Posterior of `P` (position 21) | 7.3e-5 (dense 0.99992, sparse 0.999995) |
| Root log-likelihood | 2e-13 |

With the folding threshold set to 0 (no position folds), all four cases of the test agree with dense within 1e-6, so the fold causes the difference. The case `long_branches_at_root`, the same tips with the polytomy at the root, agrees with dense, and so do both placements with branch length 0.01, where the fold error is far below 1e-6.

## Fix direction

Carry the gap and unknown references of the missing children into the candidate site table of the forward product, or keep the explicit entries of the node's message to its parent at its child-missing sites in the posterior, so that a non-root posterior keeps the same positions as the root. A position whose posterior equals the shared profile of its state still folds without loss (`fn is_shared_profile()` in `message.rs`), so the cost of the change is limited to positions that carry information.

Validation: enable the ignored case, and compare the run time of a production build with many `N` runs below polytomies, such as mpox B.1, before and after the change.

## Related issues

- [H-ancestral-sparse-dense-internal-residues-diverge.md](H-ancestral-sparse-dense-internal-residues-diverge.md): names the fold of resolved positions as one lead for the internal-node residues that differ between the representations
- [M-sparse-dense-likelihood-residual-without-missing-data.md](M-sparse-dense-likelihood-residual-without-missing-data.md): a likelihood difference that the fold may cause without missing data
