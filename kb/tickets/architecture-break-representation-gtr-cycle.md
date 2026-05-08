# Break representation/ <-> gtr/ dependency cycle

## Description

10+ files in `representation/partition/` import `GTR` from `gtr/`. Two files in `gtr/infer_gtr/` import back from `representation/`, forming a bidirectional cycle.

## Reverse direction (breakable)

- `gtr/infer_gtr/dense.rs` imports `PartitionMarginalDense` to read `(msg_to_parent, msg_to_child)` edge partition profiles
- `gtr/infer_gtr/fitch.rs` imports `PartitionCompressed` to iterate edge substitutions

## Fix

Parameterize GTR inference functions over iterators instead of concrete partition types:

- `dense.rs`: accept `impl Iterator<Item = (&Array2<f64>, &Array2<f64>)>` instead of `&PartitionMarginalDense`
- `fitch.rs`: accept `impl Iterator<Item = (usize, usize)>` (substitution count, effective length) instead of `&dyn PartitionCompressed`

Callers construct the iterator from partition data before calling GTR inference.

## Related issues

- Source: [M-core-remaining-architectural-debt-after-extraction](../issues/M-core-remaining-architectural-debt-after-extraction.md)
