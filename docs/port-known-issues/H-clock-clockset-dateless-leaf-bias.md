# ClockSet dateless-leaf contribution biases regression

`ClockSet::leaf_contribution_to_parent` at [packages/treetime/src/commands/clock/clock_set.rs#L41-L51](../../packages/treetime/src/commands/clock/clock_set.rs#L41-L51) computes weighted sums for clock regression. When `date == None` (dateless leaf, i.e. a `bad_branch` leaf), the function returns a `ClockSet` with `norm = 0` (correctly excluding the leaf from the count of dated observations), but `d_sum = branch_length / variance` and `dsq_sum = branch_length^2 / variance` are populated with non-zero values.

## Problem

During backward regression at [packages/treetime/src/commands/clock/clock_regression.rs#L85-L90](../../packages/treetime/src/commands/clock/clock_regression.rs#L85-L90), the `ClockSet` from dateless leaf edges is propagated upward and aggregated into the parent node's sum-of-children. The regression statistics that depend on the ratio of `d_sum` to `norm` are then biased:

- `clock_rate = (dt*n - t*d) / det`: immune (numerator vanishes when `t_sum = 0` for dateless leaves)
- `intercept = (d_sum - t_sum * rate) / norm`: mixes dateless `d_sum` with dated `norm`, biasing the intercept
- `chisq = 0.5 * (dsq*n - d^2 - ...) / norm`: inflated by dateless `dsq_sum`
- `r_val`: distorted for the same reason

The chi-squared distortion affects `find_best_root` ranking at [packages/treetime/src/commands/clock/find_best_root/find_best_root.rs#L64-L78](../../packages/treetime/src/commands/clock/find_best_root/find_best_root.rs#L64-L78), where the root with the lowest chi-squared is selected. Trees with many `bad_branch` leaves will have systematically inflated chi-squared values, selecting a suboptimal root.

## v0 comparison

v0 (`clock_tree.py`) uses a similar weighted-sum structure. The `bad_branch` logic is comparable: dateless leaves contribute divergence sums but not date sums. The bias exists in v0 as well, so this is not a v0/v1 divergence but rather a correctness issue in both implementations.

## Fix

For dateless leaves (`date == None`), return a zero `ClockSet` (equivalent to `outlier_contribution()`), excluding the leaf from all regression statistics. The leaf's branch length should not contribute to `d_sum`/`dsq_sum` when its date is unknown. Alternative: weight `d_sum`/`dsq_sum` by `norm` so they are zero when `norm == 0`.

## Related

- [M-clock-filter-residual-parity.md](M-clock-filter-residual-parity.md): clock filter residual differences between v0 and v1
