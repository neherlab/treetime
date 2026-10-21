# Sparse imputation of an unknown tip column ignores the posterior of the parent

With `--impute-missing-data`, sparse marginal reconstruction imputes a tip column that is unknown (`N`) from the fixed message of the most likely state of the parent, even where the parent has a posterior of its own at that column. Dense reconstruction imputes from the posterior of the parent evolved across the tip branch. The two backends therefore compute different tip posteriors at these columns, and the ambiguity report of the sparse backend misses them. [kb/decisions/ancestral-marginal-tip-reconstruction-and-imputation.md](../decisions/ancestral-marginal-tip-reconstruction-and-imputation.md) requires identical tip sequences from both backends.

## Evidence

Tree `(A:0.1,B:0.1)root;`, alignment `A` = `RA`, `B` = `NA`, JC69, imputation of the tips, ambiguity threshold 0 (states tied with the most likely state). Reconstructed tip `B`, written with the ambiguity codes of [kb/decisions/ambiguity-report-of-plausible-states.md](../decisions/ambiguity-report-of-plausible-states.md):

| Backend | `B` |
| ------- | ---- |
| dense   | `RA` |
| sparse  | `AA` |

The root posterior at column 0 ties `A` and `G`, so the posterior of `B` evolved from it also ties `A` and `G`, and the dense backend reports `R`. The sparse backend reads the fixed message of state `A`, which puts almost all mass on `A`. Both backends impute `A`, the first tied state.

## Mechanism

`fn compute_msg_to_child()` [packages/treetime/src/partition/marginal/sparse/forward.rs#L208-L226](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L208-L226) adds a variable entry to the message for the child only at columns where the child is not a non-character, so the message that a tip receives has no entry at its unknown columns. `fn leaf_posterior()` in [packages/treetime/src/partition/marginal/sparse/reconstruct.rs](../../packages/treetime/src/partition/marginal/sparse/reconstruct.rs) then falls back to the fixed message of `map_state()`, the most likely state of the parent at that column.

> [!IMPORTANT]
> **Investigation required.** Measure whether the imputed state itself differs between the backends. Under JC69 the most likely state of the evolved parent posterior equals the most likely state of the parent, so the imputed states agree in the example. Under a GTR with unequal equilibrium frequencies and a long tip branch, the evolved posterior can favor another state; run a dense/sparse comparison of imputed tips on such a case.

## Fix direction

Give the stored message of a tip edge a variable entry at every unknown column of the tip where the parent is variable: the parent posterior at that column with the message of the tip divided out, evolved across the tip branch. The sparse backend may count an unknown tip as data ([kb/issues/H-sparse-marginal-likelihood-depends-on-branch-of-unknown-leaf.md](H-sparse-marginal-likelihood-depends-on-branch-of-unknown-leaf.md)), so check what the tip contributes to the parent posterior at these columns. Keep the likelihood and the branch-length contributions unchanged, because they read the message before this step. Add the example above as a dense/sparse test of the ambiguity report.
