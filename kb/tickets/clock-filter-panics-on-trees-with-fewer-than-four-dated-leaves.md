# Clock filter panics on trees with fewer than four dated leaves

`clock_filter_inplace` at [packages/treetime/src/commands/clock/clock_filter.rs#L67-L70](../../packages/treetime/src/commands/clock/clock_filter.rs#L67-L70) computes interquartile range indices as `iq75 = (3 * n) / 4` and `iq25 = n / 4`, then unconditionally indexes into `leaf_clock_deviations[iq75]` and `leaf_clock_deviations[iq25]`. When `n == 0` (no dated leaves), the empty-slice index panics immediately. When `n < 4`, `iq75` equals the array length for some values of `n`, causing an index-out-of-bounds panic.

## Affected code

- Index computation: [packages/treetime/src/commands/clock/clock_filter.rs#L67-L70](../../packages/treetime/src/commands/clock/clock_filter.rs#L67-L70)
- Called from: [packages/treetime/src/commands/clock/run.rs#L131-L152](../../packages/treetime/src/commands/clock/run.rs#L131-L152) via `estimate_clock_model_with_prefilter`

## Trigger conditions

The function is reachable when `clock_filter_threshold > 0.0` (the default is `3.0`). Any tree where the number of leaves with valid dates is fewer than four triggers the panic. This includes:

- Trees with fewer than four leaves total
- Trees where most leaves have missing or invalid dates
- Trees where prior outlier filtering removed most dated leaves

## Reproduction

Run the clock command on a tree with fewer than four dated leaves:

```bash
treetime clock --tree=<3-leaf-tree>.nwk --dates=<dates>.tsv
```

Panics with `index out of bounds: the len is N but the index is N`.

## Fix

Guard the IQD computation with `n < 4` check. Return an error with a descriptive message ("too few dated leaves for IQD-based clock filtering; at least 4 required") rather than panicking.

## Related issues

- Source: [H-clock-filter-panic-small-trees.md](../issues/H-clock-filter-panic-small-trees.md) -- delete after full resolution
- [M-clock-filter-residual-parity.md](../issues/M-clock-filter-residual-parity.md) -- clock filter differences between v0 and v1
