# Require an optimized $T_c$ parameter for success

Make constant-$T_c$ optimization success depend on a valid optimizer result.

## Required changes

- Represent success and failure as distinct typed outcomes.
- Require a returned finite $T_c>0$ and successful termination status.
- Preserve the initial value only as input, never as an implicit successful output.
- Include optimizer cause and diagnostics in command errors.

## Validation

- Successful optimum, no best parameter, non-finite parameter, non-positive parameter, and failed termination.
- Whole-result comparisons and command-level error propagation.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/M-timetree-constant-tc-success-without-parameter.md](../issues/M-timetree-constant-tc-success-without-parameter.md)
