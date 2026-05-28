# Extract clock rate validation helper

Two timetree inference functions duplicate the same non-positive clock rate guard with identical error message.

## Current state

`timetree/inference/runner.rs:45-51` and `timetree/inference/branch_length_likelihood.rs:34-40` both check `clock_rate <= 0.0` and return `make_error!` with the same message.

## Target state

A `validate_clock_rate(rate: f64) -> Result<()>` helper in `timetree/inference/` encapsulates the check. Both sites call the helper.

## Implementation

1. Define `fn validate_clock_rate(clock_rate: f64) -> Result<()>` in `timetree/inference/` (e.g., in a shared module or one of the existing files)
2. Move the guard logic and error message into the helper
3. Replace both call sites with `validate_clock_rate(clock_rate)?;`

## Related issues

Source: `kb/issues/N-timetree-clock-rate-validation-duplication.md`
