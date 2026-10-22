# Sparse imputation leaves tips unresolved below an unknown internal node

With `--impute-missing-data`, a tip whose parent is unknown (`N`) at a column stays `N` in sparse marginal reconstruction, while dense reconstruction imputes the most likely state. This happens at every column where a whole clade is `N`, because the parent of each tip in that clade is itself unknown there ([kb/decisions/ancestral-dense-gap-rule-observed-leaf-gaps.md](../decisions/ancestral-dense-gap-rule-observed-leaf-gaps.md)). [kb/decisions/ancestral-marginal-tip-reconstruction-and-imputation.md](../decisions/ancestral-marginal-tip-reconstruction-and-imputation.md) requires that both backends produce identical tip sequences.

## Evidence

Tree `(((C1,C2)P,C3)P2,((D,E)Q,F)Q2)R` with `C1`, `C2`, `C3` = `ANGT` and `D`, `E`, `F` = `AAGT`, run with `treetime ancestral --model=jc69 --impute-missing-data --include-leaves`:

| Node        | dense  | sparse |
| ----------- | ------ | ------ |
| `P2`, `P`   | `ANGT` | `ANGT` |
| `C1`-`C3`   | `AAGT` | `ANGT` |
| other nodes | `AAGT` | `AAGT` |

The table is the output of the current code with `--representation dense` and `--representation sparse`.

## Mechanism

`fn reconstruct_leaf_sequence()` [packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L80-L105](../../packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L80-L105) calls `fn impute_state()` [packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L107-L125](../../packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L107-L125) for each unknown and ambiguous position of the tip. For a column that is fixed in the parent message, `impute_state()` looks up the fixed message by the parent's state (`fn map_state()`). At an unknown column the parent's state is `N`, the fixed message map has no `N` entry, and `impute_state()` returns the observed character, so the tip keeps its `N`. The dense backend stores a posterior profile for every column of the tip and replaces each non-canonical, non-gap character with the argmax of that profile (`fn reconstruct_leaf_sequence()` [packages/treetime/src/partition/marginal/dense/partition.rs#L208-L220](../../packages/treetime/src/partition/marginal/dense/partition.rs#L208-L220)).

## Impact

- Sparse `--impute-missing-data` output differs from dense at every column where a clade is entirely `N`
- With `--report-ambiguous`, dense reports `N5A` on the imputed tips, and sparse reports nothing

## Fix direction

Resolve the parent state before the fixed-message lookup in `impute_state()`, so that tips below an unknown parent are imputed from the message the parent passes down at that column, as the dense backend does. Add the tree above as a dense/sparse test.
