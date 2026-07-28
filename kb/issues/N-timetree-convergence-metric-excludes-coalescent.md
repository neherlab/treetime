# Convergence metric silently excludes failed coalescent log-likelihood

`TimetreeOptimizer::record()` at [packages/treetime/src/timetree/convergence/optimizer.rs#L64](../../packages/treetime/src/timetree/convergence/optimizer.rs#L64) computes `log_lh_total` by collecting `[log_lh_seq, log_lh_pos, log_lh_coal]` into an iterator, calling `.flatten()`, and reducing with addition. When `log_lh_coal` is `None` (coalescent computation failed or was skipped), `.flatten()` silently drops it from the sum.

## Impact

Trace rows where coalescent succeeded have `log_lh_total = log_lh_seq + log_lh_pos + log_lh_coal`, while rows where it failed have `log_lh_total = log_lh_seq + log_lh_pos`. The convergence delta between consecutive rows compares sums with different numbers of terms, making the convergence signal unreliable when coalescent computation is intermittent.

The per-component fields (`log_lh_seq`, `log_lh_pos`, `log_lh_coal`) are individually correct and independently comparable, so the issue is limited to the aggregate `log_lh_total` used for convergence detection.

## Affected code

- Aggregate computation: [packages/treetime/src/timetree/convergence/optimizer.rs#L64](../../packages/treetime/src/timetree/convergence/optimizer.rs#L64)
- Convergence check: `TimetreeOptimizer::next_iter()` uses `log_lh_total` delta for early stopping

## Fix

Track which components contributed to `log_lh_total` and compare only matching-component rows, or warn when the component set changes between iterations.

## Related tickets

- [kb/tickets/timetree-fix-convergence-metric-coalescent-exclusion.md](../tickets/timetree-fix-convergence-metric-coalescent-exclusion.md)
