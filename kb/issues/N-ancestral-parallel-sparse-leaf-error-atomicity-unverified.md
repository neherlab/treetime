# Parallel sparse leaf setup has no defined error atomicity contract

> [!IMPORTANT]
> **Investigation and discussion required.** This is an unverified failure-state contract, not a confirmed regression. Define the required atomicity before preparing an implementation ticket.

## Problem

`pub(crate) fn attach_seqs_to_graph()` [packages/treetime/src/ancestral/fitch.rs#L55](../../packages/treetime/src/ancestral/fitch.rs#L55) mutates leaf descriptions while resolving alignment records in parallel. If another leaf has no name or matching record, the function returns an error after an execution-dependent subset of descriptions may already have changed.

Partition construction collects a complete node map before extending each partition, which prevents partial insertion within one partition. Multiple partitions are still processed sequentially, so an error in a later partition can leave earlier partitions populated. The pre-parallel implementation also mutated state incrementally, but the parallel schedule changes which leaf-description mutations can precede an error.

Callers currently discard setup state when construction fails, which may make transactional rollback unnecessary. That ownership assumption is not stated or tested.

## Research required

- Trace every caller to establish whether a failed graph or partition set can remain observable.
- Exercise missing-name, missing-record, invalid-sequence, and multi-partition failures under multiple thread counts.
- Decide whether setup must be transactional, must leave a documented partial state, or may rely on callers discarding the state.

## Locations

- `pub(crate) fn attach_seqs_to_graph()` [packages/treetime/src/ancestral/fitch.rs#L55](../../packages/treetime/src/ancestral/fitch.rs#L55)
- Sparse marginal tests [packages/treetime/src/ancestral/__tests__/test_marginal_sparse.rs](../../packages/treetime/src/ancestral/__tests__/test_marginal_sparse.rs)

## Related issues

- [N-ancestral-parallel-sparse-leaf-validation-coverage.md](N-ancestral-parallel-sparse-leaf-validation-coverage.md)
- [N-ancestral-parallel-sparse-leaf-single-thread-regression.md](N-ancestral-parallel-sparse-leaf-single-thread-regression.md)
