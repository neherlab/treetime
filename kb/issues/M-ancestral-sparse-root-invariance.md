# Sparse root invariance violation

Felsenstein's pulley principle states that log-likelihood is invariant to root placement for reversible GTR models. [`test_prop_marginal_dense_log_lh_root_invariance`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_root_invariance_prop.rs) (`#test_prop_marginal_dense_log_lh_root_invariance`) confirms this for dense marginal reconstruction to 1e-6 (limited by `exp(Q*(t1+t2))` vs `exp(Q*t1)*exp(Q*t2)` matrix exponential path difference from edge collapse during rerooting). The sparse counterpart [`test_prop_marginal_sparse_log_lh_root_invariance`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_root_invariance_prop.rs#L44) (`#test_prop_marginal_sparse_log_lh_root_invariance`) is currently failing: sparse marginal reconstruction violates invariance by ~0.09 on the proptest minimal failing input (caterpillar tree with short branches, 4 taxa, 10 positions).

## Cause

The Fitch forward pass resolves ambiguous state sets using parent states. Under different rootings, the same unrooted tree produces different parent states at internal nodes (because the root-to-leaf direction changes). This causes:

1. Different Fitch state assignments at internal nodes
2. Different sets of edge mutations (substitutions placed on different edges)
3. Different compression patterns in the sparse representation
4. Different likelihood contributions from the marginal computation

The root parsimony SCORE is invariant (same number of mutations regardless of root), but the specific ASSIGNMENT of mutations to edges is root-dependent. Since sparse marginal operates on the Fitch-compressed representation (which encodes the specific assignment), the likelihood inherits this root dependence.

## Failing input

Proptest minimal reproducer (committed in [`proptest-regressions/commands/ancestral/__tests__/test_marginal_root_invariance_prop.txt`](../../packages/treetime/proptest-regressions/commands/ancestral/__tests__/test_marginal_root_invariance_prop.txt)):

- Tree: `(T2:0.001,(T0:0.001,(T1:0.086,T3:0.036):0.001):0.001)root:0.001;`
  (caterpillar, most branches near-zero)
- lh1=-128.716, lh2=-128.624, diff=0.092
- Rerooting target: `node_idx=66` (selects node index `66 % 2 = 0`)

Near-zero branch lengths amplify the effect: Fitch compression makes different variable/invariant classifications depending on root placement because parent states propagated through near-zero edges differ under different rootings.

## Why dense is not affected

Dense marginal operates on full probability vectors at every position, independent of Fitch compression. The only root dependence is the matrix exponential path difference from collapsing the degree-2 old root node during rerooting (`exp(Q*(t1+t2))` vs `exp(Q*t1)*exp(Q*t2)`), which is a floating-point rounding effect at ~1e-7.

## Related

- [Dense-sparse log-likelihood divergence](M-ancestral-dense-sparse-divergence.md):
  same Fitch compression mechanism causes dense and sparse to produce different likelihoods on ~2.5% of random inputs (distinct issue, same root cause).
