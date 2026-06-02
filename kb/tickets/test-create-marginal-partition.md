# Test create_marginal_partition 3-way branching

Add unit tests for `partition::create::create_marginal_partition()` which consolidates sparse, dense+infer GTR, and dense+named GTR partition creation from 4 commands.

## Test cases

- `(dense=false, model=Infer)` -> Sparse partition, inferred GTR
- `(dense=false, model=JC69)` -> Sparse partition, named GTR
- `(dense=true, model=Infer)` -> Dense partition, inferred GTR
- `(dense=true, model=JC69)` -> Dense partition, named GTR
- `(dense=None)` -> auto-detect via `infer_dense()`

Assert: partition variant, GTR alphabet, model name, root sequence consistency.

## Location

`packages/treetime/src/partition/__tests__/test_create.rs`

## Related issues

Source: `kb/issues/N-test-coverage-gaps.md`
