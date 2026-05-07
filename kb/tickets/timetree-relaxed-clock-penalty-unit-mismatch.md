# Relaxed clock unit mismatch

## Summary

The relaxed clock penalty function mixes time units (years) with evolutionary distance units (substitutions/site), inflating the penalty by approximately 110,000x on flu-rate datasets and biasing gamma values toward 1.0.

## Details

`packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs:50:`

```rust
let act_len = parent_edge.time_length().unwrap_or(opt_len)
```

`time_length()` returns years while `opt_len` is `branch_length()` in substitutions/site. The v0 reference implementation uses `clock_length` (substitutions/site) for both values, keeping the penalty dimensionally consistent.

The dimensional mismatch makes the penalty term dominate the objective for datasets with typical molecular clock rates (~0.003 subs/site/year for influenza), suppressing rate variation across branches.

## Impact

- Relaxed clock gamma estimates biased toward 1.0 (strict clock behavior)
- Rate variation across branches underestimated
- Affects all timetree runs with `--relaxed-clock` enabled

## Fix

Replace `time_length()` with `clock_length()` (= `branch_length * clock_rate`) or use a consistent unit for both `act_len` and `opt_len`.

## Related issues

- Source: [H-timetree-relaxed-clock-unit-mismatch.md](../issues/H-timetree-relaxed-clock-unit-mismatch.md) -- delete after full resolution
