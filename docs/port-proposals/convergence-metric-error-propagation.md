# Propagate convergence metric evaluation errors instead of swallowing them

When a configured coalescent model's likelihood evaluation fails, the error is downgraded to `None` and the component is excluded from `lh_total`. This makes trace rows with failed components numerically incomparable to rows where all components succeeded.

## Current behavior

At `packages/treetime/src/commands/timetree/convergence/likelihood.rs`, all three likelihood components use the same pattern:

```rust
match compute_coalescent_total_lh(graph, tc) {
  Ok(lh) => Some(lh),
  Err(e) => {
    warn!("Coalescent likelihood unavailable: {e}");
    None
  },
}
```

At `packages/treetime/src/commands/timetree/convergence/metrics.rs`:

```rust
let lh_total = [lh_seq, lh_pos, lh_coal].into_iter().flatten().reduce(|acc, v| acc + v);
```

A failed component silently disappears from `lh_total`. The tracelog CSV shows an empty cell but no indication of whether the component was disabled or failed.

## Proposed changes

### P1: Propagate errors when a configured model fails

Change `compute_coalescent_likelihood()` to return `Result<Option<f64>, Report>`:

- `Ok(None)` when no coalescent model is active (coalescent_tc is None)
- `Ok(Some(lh))` when evaluation succeeds
- `Err(e)` when a configured model fails to evaluate

`record()` propagates the error, aborting the iteration with a clear diagnostic.

### P2: Add status field to trace (alternative to P1)

Keep soft-failure behavior but add a `status` column to the tracelog CSV indicating which components contributed to `lh_total`. Downstream analysis tools can filter rows by status.

### P3: Exclude failed components from lh_total

When a configured component fails, set `lh_total = None` for that row (instead of computing a partial sum). This prevents comparing rows with different component sets.

## Context

The `warn!` upgrade (from `debug!`) in this branch makes failures visible at default log level. The `lh_total` doc comment states "sum of available components; absent components are excluded."

In practice, coalescent evaluation either works for all iterations or fails for all (the tree structure is the same each iteration). The component-set instability is a theoretical concern, not an observed failure mode.

## Related

- `docs/port-known-issues/M-timetree-coalescent-missing-leaf-and-root-contributions.md` - backward pass gap that could cause metric/inference divergence
