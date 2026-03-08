# assign_dates fails for small trees

[`assign_dates()`](../../packages/treetime/src/commands/clock/assign_dates.rs#L33)
(`#assign_dates`) has no guard when `n_leaves < MIN_GOOD_LEAVES`. The function
proceeds with insufficient data and produces unreliable results or panics.
