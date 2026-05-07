# Replace unwrap_or_default and silent fallbacks on semantically important data

## Description

Multiple locations suppress errors by substituting default values for semantically important data. A default of 0.0 for branch length can produce division-by-zero downstream or silently exclude the branch from optimization. An empty node name makes the node invisible to output serialization.

## Instances

### unwrap_or_default() on semantically important data (8 instances)

- `commands/clock/clock_filter.rs:41:` branch length defaults to 0.0
- `commands/clock/reroot.rs:178:` branch length defaults to 0.0
- `commands/clock/rtt.rs:36:` branch length defaults to 0.0
- `seq/div.rs:26:` parent divergence defaults to 0.0
- `seq/div.rs:28:` branch length defaults to 0.0
- `commands/mugration/discrete_marginal.rs:113:` node name defaults to empty string
- `commands/prune/run.rs:477-478:` indels default to empty
- `commands/timetree/optimization/relaxed_clock.rs:87:` coefficients default to zero

### branch_length().unwrap_or(one_mutation) silent fallback

`packages/treetime/src/commands/timetree/inference/runner.rs:106:`

Edges with no branch length get `one_mutation` (= 1.0 / total_sites) as a fallback when building Poisson branch-length distributions. A missing branch length could indicate a tree-loading error or an uninitialized edge, but the fallback silently assigns a plausible value.

### Composition::adjust_count saturating_add_signed silently clamps

`packages/treetime/src/seq/composition.rs:76-79:`

Uses `saturating_add_signed` which silently clamps to 0 on underflow. A negative composition count indicates a data integrity problem that should be reported, not masked.

### extract_node_times silently drops unnamed/timeless nodes

`packages/treetime/src/commands/timetree/utils.rs:63:`

Iterates all graph nodes and collects (name, time) pairs via `filter_map`. Nodes without a name or without a time are silently excluded from the convergence tracking map, making them invisible to the convergence diff count.

## Impact

- Branch length 0.0 defaults cause division-by-zero or exclusion from optimization
- Convergence tracking undercounts changes at unnamed or timeless nodes

## Related issues

- Source: [N-error-suppression-unwrap-or-defaults.md](../issues/N-error-suppression-unwrap-or-defaults.md) -- delete after full resolution
- [N-production-unwrap-expect-audit.md](../issues/N-production-unwrap-expect-audit.md) -- related unwrap/expect audit
