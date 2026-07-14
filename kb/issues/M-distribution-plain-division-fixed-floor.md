# Plain distribution division applies an unscaled fixed divisor floor

`Plain::safe_divisor()` replaces every value below `1e-10` with `1e-10`. This includes valid small positive densities, zero, and invalid negative densities. [`packages/treetime-distribution/src/policy.rs#L57-L61`](../../packages/treetime-distribution/src/policy.rs#L57-L61)

The point/function, range/function, and function/function division kernels apply this transformation before division. [`packages/treetime-distribution/src/distribution_ops/divide.rs#L55-L61`](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L55-L61) [`packages/treetime-distribution/src/distribution_ops/divide.rs#L72-L86`](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L72-L86) [`packages/treetime-distribution/src/distribution_ops/divide.rs#L108-L130`](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L108-L130)

For cavity division, let $q=2$ be the quotient to recover and $d=10^{-12}$ the divisor included in the numerator. The implemented ratio $(qd)/\max(d,10^{-10})$ returns $0.02$ instead of recovering $q=2$. The fixed floor also silently turns a negative density into a positive denominator.

## Decision required

Choose and validate the numerical contract before implementation:

- perform cavity division in negative-log space, representing zero probability explicitly; or
- reject non-finite, negative, and zero Plain divisors while dividing every strictly positive value without a fixed floor.

The decision must define zero-over-zero behavior and include a timetree cavity-message oracle. No implementation ticket is ready while this contract is unresolved.
