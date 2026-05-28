# Piecewise constant and linear function structs share identical skeleton

`PiecewiseConstantFn` (100 lines) and `PiecewiseLinearFn` (136 lines) in `treetime-grid` have identical struct layout (`breakpoints: Array1<f64>`, `values: Array1<f64>`), identical `new()` validation (ascending breakpoint check via `windows(2).all()`), identical `breakpoints()` accessor, and identical `eval_many()` tandem-sweep implementation. Only `eval()` differs (step lookup vs linear interpolation) and the value-count invariant (`n+1` constant segments vs `n` linear segments).

## Locations

- `packages/treetime-grid/src/piecewise_constant_fn.rs:1-100:` step function
- `packages/treetime-grid/src/piecewise_linear_fn.rs:1-136:` linear interpolation function

## Impact

Pure maintainability. Changes to breakpoint validation, storage, or `eval_many` require updating both files.

## Action

Extract a common `PiecewiseFnBase` struct with shared fields, validation, and `eval_many()`. `PiecewiseConstantFn` and `PiecewiseLinearFn` wrap or delegate to the base, providing only their specific `eval()` and value-count invariant.
