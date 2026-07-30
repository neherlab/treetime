# Multiplication and normalize do not compose or preserve tail metadata

`fn distribution_multiplication()` returns results with default `Error` tails regardless of operand tail policies. `fn DistributionFunction.scale_y()` (called by `fn Distribution.normalize()`) also resets tails to `Error` by constructing a new `GridFn` with defaults.

The backward pass works around this by manually re-applying `Constant` left / `Zero` right tails after each multiply+normalize step at [packages/treetime/src/timetree/inference/backward_pass.rs#L79-L82](../../packages/treetime/src/timetree/inference/backward_pass.rs#L79-L82). This is fragile: any new caller that accumulates distributions without re-applying tails will silently break under disjoint-support conditions.

## Suggested fix

Two changes:

1. `fn DistributionFunction.scale_y()` should preserve the original `left_extrap` / `right_extrap` from the source `GridFn` instead of constructing a new one with defaults. This makes `normalize()` tail-preserving.

2. `fn multiply_function_function()` (and `fn multiply_range_function()`, `fn multiply_formula_function()`) should compute result tails from operand tails. Composition rule per side:

| A tail     | B tail     | Result tail |
| ---------- | ---------- | ----------- |
| `Constant` | `Constant` | `Constant`  |
| `Constant` | `Zero`     | `Zero`      |
| `Constant` | `Error`    | `Error`     |
| `Zero`     | `Zero`     | `Zero`      |
| `Zero`     | `Error`    | `Error`     |
| `Error`    | `Error`    | `Error`     |

With both changes, the backward pass re-apply at [backward_pass.rs#L79-L82](../../packages/treetime/src/timetree/inference/backward_pass.rs#L79-L82) becomes unnecessary and can be removed. The forward pass tail application (after normalize) would still overwrite with forward-specific tails, which is correct.

## Related

- [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md): documents the current tail system
- [kb/issues/M-timetree-silent-empty-time-distribution.md](M-timetree-silent-empty-time-distribution.md): safety net for cases where tails still can't prevent empty results
