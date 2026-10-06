# Decimal-digit limit of the float formatter fails for small numbers

`fn float_to_digits()` in [packages/treetime-utils/src/fmt/float.rs](../../packages/treetime-utils/src/fmt/float.rs) passes `max_decimal_digits` to `pretty_dtoa` 0.3.0. With 3 decimal digits and no significant-digit limit:

- `-0.00001` is written as `-1.0e-5`: the limit is not applied to numbers written with an exponent
- `0.0006` panics with "attempt to subtract with overflow" at `pretty_dtoa` `lib.rs:322`. The `release` profile keeps overflow checks, so release builds panic too
- `0.00049` is written as `0`, as expected

No command sets `NwkWriteOptions::weight_decimal_digits`, so the outputs are not affected yet. The `util-newick` writer rounds to the decimal digits itself before formatting (`fn round_to_decimal_digits()` in [packages/util-newick/src/number.rs](../../packages/util-newick/src/number.rs)).

## Fix direction

Round to the decimal digits before `pretty_dtoa`, as `util-newick` does, and test values below `10^-4` and values that round up to the next digit.
