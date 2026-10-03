# Sparse imputation leaves tips unresolved below an unknown internal node

With `--impute-missing-data`, a tip whose parent is unknown (`N`) at a column stays `N` in sparse marginal reconstruction, while dense reconstruction imputes the most likely state. This happens at every column where a whole clade is `N`, because the parent of each tip in that clade is itself unknown there ([kb/decisions/ancestral-dense-gap-rule-observed-leaf-gaps.md](../decisions/ancestral-dense-gap-rule-observed-leaf-gaps.md)). [kb/decisions/ancestral-marginal-tip-reconstruction-and-imputation.md](../decisions/ancestral-marginal-tip-reconstruction-and-imputation.md) requires that both backends produce identical tip sequences.

## Evidence

Tree `(((C1,C2)P,C3)P2,((D,E)Q,F)Q2)R` with `C1`, `C2`, `C3` = `ANGT` and `D`, `E`, `F` = `AAGT`, run with `treetime ancestral --model=jc69 --impute-missing-data --include-leaves`:

| Node        | dense  | sparse |
| ----------- | ------ | ------ |
| `P2`, `P`   | `ANGT` | `ANGT` |
| `C1`-`C3`   | `AAGT` | `ANGT` |
| other nodes | `AAGT` | `AAGT` |

The `dev` branch at `3905a8f1` gives the same sparse result, so the defect predates the derived-mutation output.

## Mechanism

`fn reconstruct_leaf_sequence()` in [packages/treetime/src/partition/marginal/sparse/reconstruct.rs](../../packages/treetime/src/partition/marginal/sparse/reconstruct.rs) looks up the message from the parent for a fixed column by the parent's state (`fn map_state()`). At an unknown column the parent's state is `N`, the fixed message map has no `N` entry, and the loop skips the position, so the tip keeps its observed `N`. The dense backend evolves the parent message of every column and takes the argmax under the observed mask (`fn reconstruct_leaf_sequence()` in [packages/treetime/src/partition/marginal/dense/partition.rs](../../packages/treetime/src/partition/marginal/dense/partition.rs)).

## Impact

- Sparse `--impute-missing-data` output differs from dense at every column where a clade is entirely `N`
- With `--report-ambiguous`, dense reports `N5A` on the imputed tips, and sparse reports nothing

## Fix direction

Impute tips below an unknown parent from the message the parent passes down at that column, as the dense backend does, instead of skipping the position.
