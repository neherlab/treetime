# Approximate test assertions keep the default absolute epsilon of `approx`

`pretty_assert_ulps_eq!` [packages/treetime-utils/src/testing/assert.rs#L79-L101](../../packages/treetime-utils/src/testing/assert.rs#L79-L101) and `pretty_assert_abs_diff_eq!` [packages/treetime-utils/src/testing/assert.rs#L104-L126](../../packages/treetime-utils/src/testing/assert.rs#L104-L126) pass their options through to `approx::ulps_eq!` and `approx::abs_diff_eq!`. An option the caller leaves out takes the default of approx 0.5.1 (2022-01-23):

- **Default epsilon**: `f64::EPSILON`, about $2.2 \cdot 10^{-16}$ [[src](https://github.com/brendanzab/approx/blob/13fcbc78dcda4d0136de07c36a7979d545f08e38/src/abs_diff_eq.rs#L87)]
- **ULP comparison**: `ulps_eq` first returns true when the absolute difference is at most epsilon, and compares units in the last place only after that [[src](https://github.com/brendanzab/approx/blob/13fcbc78dcda4d0136de07c36a7979d545f08e38/src/ulps_eq.rs#L43-L47)]

So `pretty_assert_ulps_eq!(1e-20, 2e-20, max_ulps = 4)` passes although the values differ by a factor of 2, and `pretty_assert_abs_diff_eq!(a, b)` without `epsilon` checks an absolute difference of $2.2 \cdot 10^{-16}$ whatever the magnitude of the values. A check of a tiny value, such as a tail density or a small probability mass, can pass whatever the code computes.

## Fix direction

- `pretty_assert_ulps_eq!` passes `epsilon = 0.0` unless the caller gives one, so the comparison is in units in the last place only
- `pretty_assert_abs_diff_eq!` requires `epsilon`, so every absolute tolerance is visible at the call site and follows the tolerance rules of the project
- Rerun the tests after the change and treat each new failure as a finding: a wrong expected value, a tolerance to state, or a defect

## Validation

- A test of the macros: `pretty_assert_ulps_eq!(1e-20, 2e-20, max_ulps = 4)` fails
- The test suites of every crate pass, or each new failure has its own fix or issue
