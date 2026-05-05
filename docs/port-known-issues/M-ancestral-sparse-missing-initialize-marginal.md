# Sparse variant skips initialize_marginal before update_marginal

## Summary

For `--model=infer`, the sparse path calls `update_marginal` without first calling `initialize_marginal`, diverging from the dense path which calls both.

## Details

`packages/treetime/src/commands/ancestral/run.rs:105-130:`

The dense path (line 123+) calls `initialize_marginal` then `update_marginal`. The sparse path (line 105+) calls `update_marginal` only, skipping initialization. Both paths call `update_marginal` once.

`initialize_marginal` performs the first backward+forward pass with the initial GTR model and populates baseline profiles. Without it, the sparse path starts `update_marginal` from an uninitialized state, relying on whatever Fitch-era data remains in the partition.

## Impact

- Sparse marginal reconstruction with `--model=infer` may produce different results than dense
- First marginal iteration operates on different initial conditions between dense and sparse modes
- Affects ancestral reconstruction accuracy when GTR model is inferred from data

## Fix

Add `initialize_marginal` call before `update_marginal` in the sparse path, matching the dense control flow.
