# Extract common base from piecewise constant and linear function structs

Two piecewise function structs share identical struct layout, validation, and eval_many logic.

## Current state

`treetime-grid/src/piecewise_constant_fn.rs` (100 lines) and `piecewise_linear_fn.rs` (136 lines) duplicate struct fields, `new()` breakpoint validation, `breakpoints()` accessor, and `eval_many()` tandem sweep. Only `eval()` and the value-count invariant differ.

## Target state

A `PiecewiseFnBase` struct (or trait with default methods) handles shared storage, validation, and `eval_many()`. `PiecewiseConstantFn` and `PiecewiseLinearFn` provide only their specific `eval()` and value-count check.

## Implementation

1. Define `PiecewiseFnBase { breakpoints: Array1<f64>, values: Array1<f64> }` with shared validation in `new()`
2. Implement `breakpoints()`, `values()`, and `eval_many()` on the base
3. Define an `Evaluator` trait or strategy enum for the `eval()` behavior (step vs interpolation)
4. Wrap or compose the base in both concrete types, forwarding shared methods
5. Verify all existing tests pass unchanged

## Related issues

Source: `kb/issues/N-grid-piecewise-function-shared-skeleton.md`
