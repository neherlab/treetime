# Decimal-digit limit of the float formatter fails for small numbers

`fn float_to_digits()` in [packages/treetime-utils/src/fmt/float.rs](../../packages/treetime-utils/src/fmt/float.rs) passes `max_decimal_digits` to `pretty_dtoa` 0.3.0. With 3 decimal digits and no significant-digit limit:

- `-0.00001` is written as `-1.0e-5`: the limit is not applied to numbers written with an exponent
- `0.0006` panics with "attempt to subtract with overflow" at `pretty_dtoa` `lib.rs:322`. The `release` profile keeps overflow checks, so release builds panic too
- `0.00049` is written as `0`, as expected

No command sets `NwkWriteOptions::weight_decimal_digits`, so the outputs are not affected yet. The `util-newick` writer rounds to the decimal digits itself before formatting (`fn round_to_decimal_digits()` in [packages/util-newick/src/number.rs](../../packages/util-newick/src/number.rs)).

## Fix direction

Round to the decimal digits before `pretty_dtoa`, as `util-newick` does, and test values below `10^-4` and values that round up to the next digit.

## Proposed implementation plan

1. In `fn float_to_digits()` (`packages/treetime-utils/src/fmt/float.rs`), round the value to `max_decimal_digits` before `pretty_dtoa` and pass no decimal limit to `pretty_dtoa`; reuse the rounding of `fn round_to_decimal_digits()` in `packages/util-newick/src/number.rs` (move it to `treetime-utils` and call it from both, so one rounding rule owns the concern)
2. Callers that reach the path: `fn format_weight()` in `packages/treetime-io/src/nwk.rs` with `weight_decimal_digits`, and `to_decimal_digits()`/`to_digits()` of the trait in `float.rs`
3. Tests in `packages/treetime-utils/src/fmt/__tests__/` as an `rstest` table with 3 decimal digits: `-0.00001` -> `-0.0`, `0.0006` -> `0.001`, `0.00049` -> `0.0`, `0.0995` -> `0.1`, `1.2345` -> `1.235` or `1.234` by the documented tie rule. Oracle: hand calculation of decimal rounding; the `0.0006` case first fails (panic) before the fix
4. Validate with `just test-rs -E 'package(treetime-utils)'` and `-E 'package(treetime-io)'`
