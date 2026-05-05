# Clock pre-filter ignores user-provided variance parameters

## Summary

The clock command's outlier pre-filter step uses `ClockParams::default()` instead of user-provided `--sequence-length` and `--tip-slack` values, making the pre-filter insensitive to user-specified variance scaling.

## Details

`packages/treetime/src/commands/clock/run.rs:132:`

The clock command runs a pre-filter step to remove outlier tips before the main regression. This step calls `fn estimate_clock_model_with_reroot_policy` with `ClockParams::default()` instead of the user-provided parameters.

The main regression after filtering does use user parameters correctly. The pre-filter uses default variance scaling regardless of what the user specified via `--sequence-length` and `--tip-slack`.

## Impact

- Tips that should be retained (per user-specified tolerance) are incorrectly removed during pre-filtering
- Tips that should be filtered (per user-specified tolerance) pass through pre-filtering
- The effect is larger when user-specified parameters differ substantially from defaults
- Users who explicitly set `--tip-slack` expect it to govern all filtering stages

## Fix

Pass user-provided `ClockParams` (or at minimum `sequence_length` and `tip_slack`) to the pre-filter call instead of `ClockParams::default()`.
