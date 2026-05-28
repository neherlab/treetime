# ValidationRunner print methods duplicated across three runner implementations

Three `TestRunner` implementations each define 6 print methods that identically delegate to `ValidationConsole::*` static methods (18 identical method bodies total).

## Locations

- `packages/treetime-validation/src/testing/runners/multiplication.rs:46-85:` `MultiplicationRunner` - 6 methods
- `packages/treetime-validation/src/testing/runners/multiplication.rs:116-153:` `ChainMultiplicationRunner` - 6 methods
- `packages/treetime-validation/src/testing/runners/convolution.rs:44-84:` `ConvolutionRunner` - 6 methods

Delegated methods: `print_test_configuration`, `print_header`, `print_progress_table_header`, `print_success_row`, `print_failure_row`, `print_error_summary`.

## Impact

Pure maintainability. Adding a runner requires copying 6 trivial delegation methods.

## Action

Add default method implementations to the `TestRunner` trait for all 6 print methods, each delegating to the corresponding `ValidationConsole::*` method.
