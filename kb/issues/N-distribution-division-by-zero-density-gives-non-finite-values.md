# Division by a divisor with zero density gives non-finite values

`fn distribution_division()` [packages/treetime-distribution/src/distribution_ops/divide.rs#L19-L42](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L19-L42) subtracts neg-log ordinates (`NegLog::divide` is $a - b$, [packages/treetime-distribution/src/policy.rs#L91-L93](../../packages/treetime-distribution/src/policy.rs#L91-L93)), and `NegLog::safe_divisor` passes the divisor through unchanged ([packages/treetime-distribution/src/policy.rs#L99-L101](../../packages/treetime-distribution/src/policy.rs#L99-L101)). Where the divisor has zero density ($+\infty$) at an evaluated point, three branches produce values that are not a density:

- **Range by function** (`fn divide_range_by_function()` [packages/treetime-distribution/src/distribution_ops/divide.rs#L79-L108](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L79-L108)): the finite range amplitude minus $+\infty$ gives $-\infty$, an infinite density, at every grid point where the divisor is zero; the branch checks no value
- **Single-point intersection** of a range or a function with a function (`SupportIntersection::Point` in both functions): the dividend value minus $+\infty$ gives $-\infty$, and $+\infty - (+\infty)$ gives `NaN` when the dividend is zero there too
- **Point by function** (`fn divide_point_by_function()` [packages/treetime-distribution/src/distribution_ops/divide.rs#L44-L77](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L44-L77)) checks its value and returns the guarded empty result, so it does not produce a non-finite value

`fn divide_function_by_function()` handles a shared zero-density end point of a function dividend and divisor through `fn shared_hard_end()`, which compares the power laws of both ends. The branches above have no such rule.

The forward pass divides only function posteriors by function messages ([packages/treetime/src/timetree/inference/forward_pass.rs#L308](../../packages/treetime/src/timetree/inference/forward_pass.rs#L308)), so no measured run reaches these branches.

## Fix direction

Give each branch the rule of `fn shared_hard_end()` or an error that names the point: a zero-density divisor under a nonzero dividend is a division by zero, and a shared zero is the limit of the two laws.

## Validation

Tests through `fn distribution_division()`: a range divided by a function with a zero-density grid point inside the range, and a single-point intersection where the divisor is zero, each either return a defined result or the error.
