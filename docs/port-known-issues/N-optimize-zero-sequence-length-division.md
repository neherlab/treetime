# Division by zero when total sequence length is zero

`run_optimize_mixed()` and `initial_guess_mixed()` compute `one_mutation = 1.0 / total_length as f64` where `total_length` is the sum of sequence lengths across all partitions. If all partitions have zero sequence length, this produces `+inf`, which propagates into grid search bounds, Newton tolerance, and initial guess values.

## Location

- [packages/treetime/src/commands/optimize/optimize_unified.rs#L357](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L357): `let one_mutation = 1.0 / total_length as f64;` in `run_optimize_mixed()`
- [packages/treetime/src/commands/optimize/optimize_unified.rs#L473](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L473): same in `initial_guess_mixed()`

## Impact

Medium. An empty alignment (zero sequence length) is invalid input that should be rejected upstream. No current code path reaches this with zero-length partitions. The risk is future code paths or malformed input bypassing validation.

## Fix

Add an early check: `if total_length == 0 { return make_error!("...") }` before computing `one_mutation`.
