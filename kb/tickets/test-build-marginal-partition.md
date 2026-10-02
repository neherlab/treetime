# Test build_marginal_partition for each representation and GTR model

Add unit tests for `partition::create::build_marginal_partition()` which consolidates sparse, dense+infer GTR, and dense+named GTR partition creation for the ancestral, optimize and timetree commands (prune calls `build_sparse_partition()` directly).

## Test cases

- `(Representation::Sparse, model=Infer)` -> Sparse partition, inferred GTR
- `(Representation::Sparse, model=JC69)` -> Sparse partition, named GTR
- `(Representation::Dense, model=Infer)` -> Dense partition, inferred GTR
- `(Representation::Dense, model=JC69)` -> Dense partition, named GTR

Assert: partition variant, GTR alphabet, root sequence consistency. The `--dense` default (`Representation::resolve(None)` through `infer_dense()`) is covered by the plan-resolution tests in `packages/treetime/src/ancestral/__tests__/test_plan.rs`.

## Location

`packages/treetime/src/partition/__tests__/test_create.rs`

## Related issues

Source: [kb/issues/N-test-coverage-gaps.md](../issues/N-test-coverage-gaps.md)
