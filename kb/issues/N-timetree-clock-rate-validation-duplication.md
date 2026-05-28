# Clock rate non-positive guard duplicated at two sites

Identical `if clock_rate <= 0.0` guard with the same error message appears in two timetree inference functions. Both return `make_error!` with the same multi-line message about negative clock rate and suggesting `--clock-rate`.

## Locations

- `packages/treetime/src/timetree/inference/runner.rs:45-51:` guard in `run_timetree_inference()`
- `packages/treetime/src/timetree/inference/branch_length_likelihood.rs:34-40:` guard in branch length likelihood computation

Both have identical code:

```rust
if clock_rate <= 0.0 {
  return make_error!(
    "Clock rate estimate is negative ({clock_rate:.6e}). \
     The data may lack sufficient temporal signal for automatic rate estimation. \
     Please specify --clock-rate explicitly."
  );
}
```

## Impact

Trivial duplication. If the validation logic or message changes, two sites need updating.

## Action

Extract `validate_clock_rate(rate: f64) -> Result<()>` helper in `timetree/inference/`.
