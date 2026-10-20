# Distribution Range accepts invalid supports and amplitudes

> [!IMPORTANT]
> **Decision required.** The range construction contract is unapproved: fallible range construction with validated endpoints and amplitude, or a validated support type that separates positive-width intervals from points. See "Decision required" below.

`Distribution::range()` is infallible and delegates to `DistributionRange::new()`. [`packages/treetime-distribution/src/distribution_core/distribution.rs#L42-L44`](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L42-L44)

`DistributionRange::new()` stores arbitrary endpoints and amplitude without checking ordering, finiteness, or the `YAxisPolicy` validity rule. [`packages/treetime-distribution/src/distribution_core/range.rs#L15-L21`](../../packages/treetime-distribution/src/distribution_core/range.rs#L15-L21)

Reversed supports can flow into operations as if they were empty or partially overlapping. NaN endpoints are more dangerous: `fn f64::min()` and `fn f64::max()` return the non-NaN operand when exactly one operand is NaN [[doc](https://doc.rust-lang.org/std/primitive.f64.html#method.min)], turning an invalid support into a valid-looking intersection and concealing the original numerical failure.

The production callers pass ordered, finite endpoints (date ranges, negation, overlaps, convolution of ranges), so no runtime failure has been shown.

## Decision required

Choose one construction contract:

- make range construction fallible and validate finite ordered endpoints plus policy-valid amplitude; or
- introduce a validated support type that distinguishes positive-width intervals from points before constructing a distribution.

Update every direct constructor consistently. Implementation waits until this API decision is approved.
