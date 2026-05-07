# Polytomy resolution leaves stale fields after topology change

## Summary

`prepare_tree_after_topology_change` clears only `msg_to_parent` and `branch_length_distribution` after polytomy resolution, leaving `gamma`, `clock_to_parent`, `clock_to_child`, and `clock_from_child` stale on newly created edges.

## Details

`packages/treetime/src/commands/timetree/optimization/polytomy.rs:456-478:`

After resolving a polytomy (splitting a multifurcation into binary nodes), the function resets the message-passing fields but leaves relaxed-clock and time-tree specific fields at their values from the previous iteration:

- `gamma`: relaxed clock rate multiplier (should be 1.0 for new edges)
- `clock_to_parent`: clock-scaled branch length toward parent
- `clock_to_child`: clock-scaled branch length toward child
- `clock_from_child`: clock message from child

These stale values from the pre-resolution topology persist into the next refinement iteration, biasing the optimizer with rate information from edges that no longer exist in the same configuration.

## Impact

- Relaxed clock estimates on newly resolved edges start from stale gamma values
- Time inference on resolved polytomies uses outdated clock messages
- Effects compound across refinement iterations

## Fix

Reset `gamma` to 1.0 and clear `clock_to_parent`, `clock_to_child`, `clock_from_child` alongside the existing field resets in `prepare_tree_after_topology_change`.

## Related issues

- Source: [M-timetree-polytomy-stale-fields-after-topology-change.md](../issues/M-timetree-polytomy-stale-fields-after-topology-change.md) -- delete after full resolution
