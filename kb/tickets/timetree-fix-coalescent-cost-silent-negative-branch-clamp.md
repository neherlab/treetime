# Fix sum_coalescent_cost silent negative branch length clamp

`sum_coalescent_cost` at [packages/treetime/src/commands/timetree/coalescent/edge_data.rs#L121](../../packages/treetime/src/commands/timetree/coalescent/edge_data.rs#L121) uses `branch_length.max(0.0)` to clamp negative branch lengths to zero. When a parent node has a younger time than its child (a transient state during optimization), the clamped branch length produces `t_merger == t_node`, so the integral `I(t_merger) - I(t_node) = 0` and the `-log(lambda)` merger cost is still added.

## Impact

This silent clamping hides parent-younger-than-child topology-time inconsistencies from the coalescent objective function. The optimizer sees a flat cost for all negative branch lengths, losing the gradient signal that would push the solution toward consistent node times. The coalescent cost appears finite and well-behaved even when the tree topology is temporally inverted.

## Affected code

- Clamp site: [packages/treetime/src/commands/timetree/coalescent/edge_data.rs#L121](../../packages/treetime/src/commands/timetree/coalescent/edge_data.rs#L121)

## Fix

Replace the silent clamp with either:

- An explicit error return when `branch_length < 0` (strict: forces callers to handle temporal inconsistency)
- A large penalty term added to the cost (soft: preserves optimizer continuity while penalizing inconsistency)
- A warning log with the current clamp behavior preserved (minimal: documents the situation without changing behavior)

## Related issues

- Source: [M-timetree-coalescent-branch-length-clamp.md](../issues/M-timetree-coalescent-branch-length-clamp.md) -- delete after full resolution
- [N-core-branch-length-clamping.md](../issues/N-core-branch-length-clamping.md) -- branch-length clamping hack in marginal passes (different code path, same pattern)
- [M-optimize-negative-branch-length-validation.md](../issues/M-optimize-negative-branch-length-validation.md) -- missing negative branch-length validation in the optimizer
