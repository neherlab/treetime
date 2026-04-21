# Clock covariation divides by zero when sequence length is absent

The `--covariation` branch in `run_clock` at [packages/treetime/src/commands/clock/run.rs#L61-L69](../../packages/treetime/src/commands/clock/run.rs#L61-L69) reads `seq_len = sequence_length.unwrap_or(0) as f64` and then uses it as a denominator in `variance_factor = overdispersion / seq_len` and `variance_offset_leaf = tip_slack^2 / seq_len / seq_len`. With the default `None` (no `--sequence-length` flag), `seq_len = 0.0` and both fields become `inf`.

## Affected code

- Division: [packages/treetime/src/commands/clock/run.rs#L61-L69](../../packages/treetime/src/commands/clock/run.rs#L61-L69)
- The inline comment on line 62 acknowledges the defect: "should error if None and covariation is true"

## Impact

Downstream `ClockParams` fields `variance_factor` and `variance_offset_leaf` are `inf`, causing:

- All regression variances to become `inf`
- `ClockSet` weighted sums to be zero (all terms divided by `inf`)
- Clock model output to be degenerate (rate = 0, intercept = NaN or 0)
- `find_best_root` chi-squared comparisons to be meaningless

The command does not crash but produces silently wrong results.

## Fix

Add a `make_error!` guard before the division: if `covariation` is true and `sequence_length` is `None`, return an error requiring the user to provide `--sequence-length`.
