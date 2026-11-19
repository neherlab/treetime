# Plain distribution division applies an unscaled fixed divisor floor

> [!IMPORTANT]
> **Decision required.** The numerical contract of Plain division is unresolved: cavity division in negative-log space, or rejection of non-finite, negative, and zero divisors without a fixed floor. See "Decision required" below.

`Plain::safe_divisor()` replaces every value below `1e-10` with `1e-10`. This includes valid small positive densities, zero, and invalid negative densities. [`packages/treetime-grid/src/y_axis_policy.rs#L49-L52`](../../packages/treetime-grid/src/y_axis_policy.rs#L49-L52)

The point/function, range/function, and function/function division kernels apply this transformation before division. [`packages/treetime-distribution/src/distribution_ops/divide.rs#L66-L67`](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L66-L67) [`packages/treetime-distribution/src/distribution_ops/divide.rs#L92-L101`](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L92-L101) [`packages/treetime-distribution/src/distribution_ops/divide.rs#L123-L134`](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L123-L134)

For cavity division, let $q=2$ be the quotient to recover and $d=10^{-12}$ the divisor included in the numerator. The implemented ratio $(qd)/\max(d,10^{-10})$ returns $0.02$ instead of recovering $q=2$. The fixed floor also silently turns a negative density into a positive denominator.

A zero divisor density at an end of the function/function quotient grid does not reach the floor: `fn shared_hard_end()` replaces that end value under both policies by the limit of the two power laws that `fn GridFn::hard_cell()` fits there (for example $[0, 1, 4] / [0, 1, 4]$ on $[0, 2]$ stores $1$ at $t = 0$, the limit of $t^2 / t^2$). Zeros inside the grid, and the point/function and range/function kernels, still use the floor.

No production path reaches the floor. The only production division is the cavity division of the timetree forward pass [`packages/treetime/src/timetree/inference/forward_pass.rs#L308`](../../packages/treetime/src/timetree/inference/forward_pass.rs#L308), which runs on `Distribution<NegLog>`. `NegLog::safe_divisor()` returns its argument unchanged [`packages/treetime-grid/src/y_axis_policy.rs#L98-L100`](../../packages/treetime-grid/src/y_axis_policy.rs#L98-L100). The floor applies only to callers of `distribution_division()` with the `Plain` policy.

## Decision required

Choose and validate the numerical contract before implementation:

- perform cavity division in negative-log space, representing zero probability explicitly; or
- reject non-finite, negative, and zero Plain divisors while dividing every strictly positive value without a fixed floor.

The decision must define zero-over-zero behavior and include a timetree cavity-message oracle. Implementation waits until this contract is resolved.
