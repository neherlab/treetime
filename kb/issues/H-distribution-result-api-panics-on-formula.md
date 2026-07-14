# Result-returning distribution operations panic on Formula operands

Several public distribution operations return `Result<Distribution<_>, Report>` but unconditionally panic when a valid `Distribution::Formula` reaches an unsupported dispatch arm:

- division: [`packages/treetime-distribution/src/distribution_ops/divide.rs#L9-L50`](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L9-L50);
- mapping: [`packages/treetime-distribution/src/distribution_ops/map.rs#L5-L19`](../../packages/treetime-distribution/src/distribution_ops/map.rs#L5-L19);
- convolution: [`packages/treetime-distribution/src/distribution_ops/convolve.rs#L13-L42`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L13-L42);
- scalar multiplication: [`packages/treetime-distribution/src/distribution_ops/scalar_multiply.rs#L9-L32`](../../packages/treetime-distribution/src/distribution_ops/scalar_multiply.rs#L9-L32).

`Formula` is part of the public `Distribution` enum and is produced by coalescent calculations. The type signatures neither exclude formulas nor state an unwind precondition. A valid enum value can therefore terminate the process instead of returning an actionable error.

## Required behavior

Every `Result`-returning distribution operation must return `Ok` or `Err` for every public distribution variant. Unsupported formula combinations should return a contextual error until their mathematical implementation is specified. Adding analytic implementations is an independent capability and must preserve the approved support and zero-divisor semantics.

## Related tickets

- [kb/tickets/safety-return-errors-for-unsupported-formula-operations.md](../tickets/safety-return-errors-for-unsupported-formula-operations.md)
