# Parallel sparse leaf setup has no defined error atomicity contract

> [!IMPORTANT]
> **Investigation and discussion required.** This is an unverified failure-state contract, not a confirmed regression. Define the required atomicity before preparing an implementation ticket.

## Problem

`pub(crate) fn attach_seqs_to_graph()` mutates each leaf description inside a fallible parallel iterator before collecting the complete result. If another leaf has no name or matching alignment record, descriptions written by successful workers remain observable even though the function returns an error.

Sequence conversion is staged into a complete node map before one partition is extended, so a conversion error cannot partially populate that partition. Partitions are processed sequentially, however: conversion failure in a later partition leaves earlier partitions populated. Edge payloads are inserted only after every partition's nodes succeed.

The function borrows the graph and partitions rather than consuming a disposable builder, so callers can retain the mutated values after failure. No API contract or test establishes that callers always discard them.

## Decision axes

### A1. Atomicity scope

- O1. Make the complete call transactional across leaf descriptions, every partition's nodes, and edge payloads. Any error leaves all inputs unchanged.
- O2. Guarantee atomicity per partition while documenting that earlier partitions and leaf descriptions may be committed. This exposes call-order-dependent state and requires every caller to reason about partial initialization.
- O3. Move setup behind a consuming builder whose error path cannot return the graph or partitions. Partial internal mutation becomes unobservable, but this changes the construction boundary and still requires proof that no shared handles escape.

**Recommendation:** O1. Stage validated leaf metadata and every partition's node and edge maps without mutating shared state, then commit all staged values after all fallible work succeeds.

### A2. Parallel failure behavior

- O1. Preserve parallel validation and construction, collect results in deterministic graph and partition order, and perform one serial commit.
- O2. Validate serially before parallel construction. This avoids concurrent side effects but duplicates traversal and does not by itself make multi-partition commit atomic.

**Recommendation:** O1. Parallel workers should produce owned staged values only; shared-state mutation belongs in the infallible commit step.

## Locations

- `pub(crate) fn attach_seqs_to_graph()` mutates descriptions during collection [packages/treetime/src/ancestral/fitch.rs#L55-L91](../../packages/treetime/src/ancestral/fitch.rs#L55-L91)
- Partition node and edge commits [packages/treetime/src/ancestral/fitch.rs#L93-L109](../../packages/treetime/src/ancestral/fitch.rs#L93-L109)
- Sparse marginal tests [packages/treetime/src/ancestral/__tests__/test_marginal_sparse.rs](../../packages/treetime/src/ancestral/__tests__/test_marginal_sparse.rs)

## Validation

- Inject missing-name, missing-record, and invalid-sequence failures at each leaf position under one and multiple worker threads.
- Use at least two partitions and force the second partition to fail conversion; compare the complete graph and every partition with their pre-call values.
- Verify successful setup produces identical deterministic state across worker counts.

## Related issues

- [N-ancestral-parallel-sparse-leaf-validation-coverage.md](N-ancestral-parallel-sparse-leaf-validation-coverage.md)
- [N-ancestral-parallel-sparse-leaf-single-thread-regression.md](N-ancestral-parallel-sparse-leaf-single-thread-regression.md)
- [M-inference-fallible-parallel-passes-partially-commit.md](M-inference-fallible-parallel-passes-partially-commit.md)
