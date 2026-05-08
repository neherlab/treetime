# Remove dead PartitionMarginal marker trait

## Description

`PartitionMarginal` in `representation/partition/traits.rs:115:` is an empty marker trait. Implemented as empty `impl` on `PartitionMarginalDense` and `PartitionMarginalSparse`. Used only as a supertrait bound on `PartitionMarginalOps`. Adds no methods, no type discrimination.

## Fix

Remove the trait definition. Remove both empty impls. Drop the supertrait bound from `PartitionMarginalOps`.

3 deletions, 1 edit.

## Related issues

- Source: [M-core-remaining-architectural-debt-after-extraction](../issues/M-core-remaining-architectural-debt-after-extraction.md)
