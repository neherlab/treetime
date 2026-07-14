# Fix convergence metric silently excluding failed coalescent likelihood

`TimetreeOptimizer::record()` at [packages/treetime/src/timetree/convergence/optimizer.rs#L64](../../packages/treetime/src/timetree/convergence/optimizer.rs#L64) computes `lh_total` by collecting `[lh_seq, lh_pos, lh_coal]` into an iterator, calling `.flatten()`, and reducing with addition. When `lh_coal` is `None` (coalescent computation failed or was skipped), `.flatten()` silently drops it from the sum.

## Impact

Trace rows where coalescent succeeded have `lh_total = lh_seq + lh_pos + lh_coal`, while rows where it failed have `lh_total = lh_seq + lh_pos`. The convergence delta between consecutive rows compares sums with different numbers of terms, making the convergence signal unreliable when coalescent computation is intermittent.

The per-component fields (`lh_seq`, `lh_pos`, `lh_coal`) are individually correct and independently comparable, so the issue is limited to the aggregate `lh_total` used for convergence detection.

## Affected code

- Aggregate computation: [packages/treetime/src/timetree/convergence/optimizer.rs#L64](../../packages/treetime/src/timetree/convergence/optimizer.rs#L64)
- Convergence check: `TimetreeOptimizer::next_iter()` uses `lh_total` delta for early stopping

## Fix

Track which components contributed to `lh_total` and compare only matching-component rows, or warn when the component set changes between iterations.

## Related issues

- Source: [kb/issues/N-timetree-convergence-metric-excludes-coalescent.md](../issues/N-timetree-convergence-metric-excludes-coalescent.md) -- delete after full resolution
