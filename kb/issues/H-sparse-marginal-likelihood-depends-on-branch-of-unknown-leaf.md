# Sparse marginal likelihood depends on the branch length of an unknown leaf

## Summary

A leaf whose sequence is unknown at every position (all `N`) carries no data, so the tree likelihood must not depend on the length of its branch. The dense marginal likelihood is independent of it, as expected. The sparse marginal likelihood changes with that branch length and differs from the dense value at the same branch lengths.

## Reproduction

Tree `((A:0.1,X:x)n1:0.1,(B:0.1,Y:0.1)n2:0.1)root;` with JC69 and the alignment

```
>A
ACGTACGTACGTACGTACGT
>B
ACGTACGTACGTACGTACAA
>X
NNNNNNNNNNNNNNNNNNNN
>Y
NNNNNNNNNNNNNNNNNNNN
```

Marginal log-likelihood after `marginal_update()` with the setup functions `setup_dense_reconstruction()` and `setup_reconstruction()` of `packages/treetime/src/optimize/__tests__/test_convergence/test_convergence_support.rs`:

| Backend | $x = 0.1$ | $x = 2.0$ |
| ------- | --------- | --------- |
| Dense   | -38.9449  | -38.9449  |
| Sparse  | -44.6571  | -65.8681  |

The same values appear on `rust` at `3825359a2` with `treetime optimize --model=jc69 --max-iter=1` on this tree. In the second iteration the branch lengths of the two backends have moved apart, and the log-likelihoods at $x = 0.1$ are -36.433 (sparse) and -48.346 (dense).

The sparse value falls with $x$, so the sparse backend treats the unknown leaf as if it had a known sequence.

> [!IMPORTANT]
> **Investigation required.** The cause is not traced. Candidates: the Fitch compression or the fixed-position composition of the sparse partition counts the unknown leaf at positions without variation, or the sparse leaf messages do not mark unknown ranges as missing data. Trace a leaf with all positions unknown through `create_fitch_partition()` and the sparse marginal messages, and check whether partly unknown leaves are affected the same way.

## Consequences

- Sparse branch lengths and likelihoods are biased for datasets with leaves that have long unknown ranges
- The default path produces such leaves: `fn complete_alignment_for_leaves()` fills every tip without a sequence with unknown characters and only warns, as long as no more than a third of the tips lack a sequence [packages/treetime/src/ancestral/attach.rs#L15-L53](../../packages/treetime/src/ancestral/attach.rs#L15-L53). Sparse is the default representation, so default runs of every sequence command are affected whenever a tip has no sequence
- `test_run_optimize_loop_returns_first_iteration_when_likelihood_declines_from_start` in `packages/treetime/src/optimize/__tests__/test_run_optimize_loop.rs` uses leaves without data and therefore runs on the dense backend only

## Validation

- A test asserts that the sparse marginal likelihood of the tree above does not depend on $x$ and equals the dense value
- The loop test above runs on the sparse backend too

## Related issues

The sparse backend is the default of every run ([N-representation-infer-dense-stub.md](N-representation-infer-dense-stub.md)), and the dense backend is its independent oracle. These issues show the two backends disagreeing on default settings and may share one cause. The tree with one unknown leaf of this issue is the smallest reproducing case, so an investigation of the others starts here.

- [H-timetree-sparse-dense-node-times-differ.md](H-timetree-sparse-dense-node-times-differ.md)
- [H-ancestral-sparse-dense-internal-residues-diverge.md](H-ancestral-sparse-dense-internal-residues-diverge.md)
- [M-ancestral-sparse-root-invariance.md](M-ancestral-sparse-root-invariance.md)
