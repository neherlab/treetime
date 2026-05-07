# Replace 19+ element-by-element float comparison loops with array assertions

## Problem

`packages/treetime/src/alphabet/__tests__/` -- 19+ instances of `for (actual, expected) in ... { assert_ulps_eq!(...) }`. Should use `pretty_assert_ulps_eq!` or `pretty_assert_abs_diff_eq!` on whole arrays.

## Related issues

- Source: [N-test-quality-deficiencies.md](../issues/N-test-quality-deficiencies.md) -- delete after full resolution
