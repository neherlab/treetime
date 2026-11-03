# FFT convolver computes the transform length without overflow checks

`FftConvolver::new()` at [packages/treetime-ops/src/convolution.rs#L87-L89](../../packages/treetime-ops/src/convolution.rs#L87-L89) checks only the minimum operand lengths, then computes `f_len + g_len - 1`. `smooth_even_len()` at [packages/treetime-ops/src/convolution.rs#L248-L258](../../packages/treetime-ops/src/convolution.rs#L248-L258) multiplies `2 * odd * ...` and calls `next_power_of_two()` without checks. For lengths near `usize::MAX`, a debug build panics and a release build wraps to a short transform.

The callers pass the lengths of allocated arrays, which the grid point limit of `MaxGridPoints` bounds far below these values, so no input of TreeTime reaches the overflow. The function is public, and a caller with other lengths gets a panic instead of an error.

## Fix

Compute the output length and the transform length with checked arithmetic, and return an error that names both operand lengths when the transform length cannot be represented.

## Validation

A test of `FftConvolver::new(usize::MAX, 2, &FftPlans::new())` asserts the error message with `assert_error!`.
