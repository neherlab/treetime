# Add default methods to TestRunner trait for print delegation

Three TestRunner implementations each define 6 identical print methods delegating to ValidationConsole.

## Current state

`MultiplicationRunner`, `ChainMultiplicationRunner`, and `ConvolutionRunner` each implement `print_test_configuration`, `print_header`, `print_progress_table_header`, `print_success_row`, `print_failure_row`, `print_error_summary` with identical bodies.

## Target state

The `TestRunner` trait provides default implementations for all 6 print methods. Runner implementations inherit them without boilerplate.

## Implementation

1. Add default method bodies to the `TestRunner` trait definition, each delegating to the corresponding `ValidationConsole::*` method
2. Remove the 18 explicit method implementations from the three runner structs
3. Verify test output unchanged

## Related issues

Source: `kb/issues/N-validation-runner-print-method-boilerplate.md`
