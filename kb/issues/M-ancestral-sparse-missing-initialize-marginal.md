# Missing initialize_marginal before update_marginal in sparse and AA dense paths

## Summary

Two code paths call `update_marginal` without first calling `initialize_marginal`, diverging from the nuc dense path which calls both.

## Affected paths

### Nuc sparse path

`packages/treetime/src/commands/ancestral/run.rs` (pipeline.rs sparse branch):

The nuc dense path calls `initialize_marginal` then `update_marginal`. The nuc sparse path calls `update_marginal` only, skipping initialization.

### AA multi-partition path (dense with named model)

`packages/treetime/src/ancestral/multi.rs` (`reconstruct_marginal_partitions`):

Calls `create_marginal_partition` then `attach_sequences` then `update_marginal`, but never `initialize_marginal`. For `--aa-model infer` (default) this is safe because Fitch bootstrapping inside `create_marginal_partition` populates initial profiles. For `--aa-model jtt92 --dense=true`, the partition is created via `PartitionMarginalDense::new()` with empty profiles, so `initialize_marginal` is needed to seed baseline marginal profiles before `update_marginal`.

## Details

`initialize_marginal` performs the first backward+forward pass with the initial GTR model and populates baseline profiles. Without it, `update_marginal` starts from an uninitialized state, relying on whatever Fitch-era data remains in the partition (nuc sparse) or empty arrays (AA dense with named model).

## Impact

- Sparse marginal reconstruction with `--model=infer` may produce different results than dense
- AA reconstruction with `--aa-model jtt92 --dense=true` operates on empty initial profiles
- First marginal iteration operates on different initial conditions between affected and unaffected modes

## Fix

Add `initialize_marginal` call before `update_marginal` in both the nuc sparse path and the AA multi-partition path, matching the nuc dense control flow.
