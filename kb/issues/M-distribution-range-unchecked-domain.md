# Distribution Range accepts invalid supports and amplitudes

`Distribution::range()` is infallible and delegates to `DistributionRange::new()`. [`packages/treetime-distribution/src/distribution_core/distribution.rs#L36-L42`](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L36-L42)

`DistributionRange::new()` stores arbitrary endpoints and amplitude without checking ordering, finiteness, or the `YAxisPolicy` validity rule. [`packages/treetime-distribution/src/distribution_core/range.rs#L17-L24`](../../packages/treetime-distribution/src/distribution_core/range.rs#L17-L24)

Reversed supports can flow into operations as if they were empty or partially overlapping. NaN endpoints are more dangerous: `fn f64::min()` and `fn f64::max()` return the non-NaN operand when exactly one operand is NaN [[doc](https://doc.rust-lang.org/std/primitive.f64.html#method.min)], turning an invalid support into a valid-looking intersection and concealing the original numerical failure.

## Decision required

Choose one construction contract:

- make range construction fallible and validate finite ordered endpoints plus policy-valid amplitude; or
- introduce a validated support type that distinguishes positive-width intervals from points before constructing a distribution.

Update serde deserialization and every direct constructor consistently. No implementation ticket is ready until this API decision is approved.
