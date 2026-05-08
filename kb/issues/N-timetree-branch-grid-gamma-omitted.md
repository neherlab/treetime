# Branch grid extent uses base clock_rate, not effective rate

`create_simple_grid()` computes the grid extent as `MAX_BRANCH_TIME * clock_rate`, but the time conversion in `compute_branch_length_distribution()` uses `effective_clock_rate = clock_rate * gamma`. The actual time coverage is `MAX_BRANCH_TIME / gamma`, not `MAX_BRANCH_TIME`.

For gamma = 2.0, coverage is 100 years. For gamma = 0.5, coverage is 400 years. Both are adequate for practical datasets. The concern is theoretical for extreme gamma values (>10), which would reduce coverage below 20 years.

v0 has the same pattern: `MAX_BRANCH_LENGTH = 4.0` is fixed in subs/site space, and gamma is applied separately to the time conversion.

## Fix

Pass gamma to `create_simple_grid` and compute `rate_max_bl = MAX_BRANCH_TIME * clock_rate * gamma`.
Or document the `MAX_BRANCH_TIME / gamma` coverage in the `MAX_BRANCH_TIME` doc comment.
