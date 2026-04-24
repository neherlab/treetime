# Shared graph topology operations misplaced inside command modules

Graph topology operations used by multiple commands are defined inside `commands/clock/` rather than in a shared module. This creates cross-command coupling, independent reimplementations, and architectural inconsistency with the existing `topology_cleanup/` extraction effort.

## Problem

[packages/treetime/src/commands/clock/reroot.rs](../../packages/treetime/src/commands/clock/reroot.rs) (288 lines) contains generic graph operations that are not clock-specific:

- `reroot_in_place()` [packages/treetime/src/commands/clock/reroot.rs#L88-L154](../../packages/treetime/src/commands/clock/reroot.rs#L88-L154) -- orchestrates graph rerooting: find best root, split edge, invert path, remove trivial old root
- `create_new_root_node()` [packages/treetime/src/commands/clock/reroot.rs#L158-L199](../../packages/treetime/src/commands/clock/reroot.rs#L158-L199) -- split an edge by inserting a new node at a fractional position
- `apply_reroot()` [packages/treetime/src/commands/clock/reroot.rs#L202-L234](../../packages/treetime/src/commands/clock/reroot.rs#L202-L234) -- invert edges along the path from old root to new root
- `remove_node_if_trivial()` [packages/treetime/src/commands/clock/reroot.rs#L239-L288](../../packages/treetime/src/commands/clock/reroot.rs#L239-L288) -- remove a degree-2 node, merging its parent and child edges

Five data types are defined alongside these functions:

- `RerootParams` [packages/treetime/src/commands/clock/reroot.rs#L17-L33](../../packages/treetime/src/commands/clock/reroot.rs#L17-L33)
- `EdgeSplitInfo` [packages/treetime/src/commands/clock/reroot.rs#L36-L48](../../packages/treetime/src/commands/clock/reroot.rs#L36-L48)
- `EdgeMergeInfo` [packages/treetime/src/commands/clock/reroot.rs#L51-L61](../../packages/treetime/src/commands/clock/reroot.rs#L51-L61)
- `RerootResult` [packages/treetime/src/commands/clock/reroot.rs#L64-L72](../../packages/treetime/src/commands/clock/reroot.rs#L64-L72)
- `RerootChanges` [packages/treetime/src/commands/clock/reroot.rs#L77-L86](../../packages/treetime/src/commands/clock/reroot.rs#L77-L86)

Only `reroot_in_place()` has a clock-specific dependency: it calls `find_best_root()`, which uses clock regression to score root positions. The remaining functions and all five data types operate on generic `Graph<N, E, D>` with `GraphNode + GraphEdge` bounds and contain no clock logic.

## Cross-command coupling

The `timetree` command imports 48 items from `commands/clock/`, many for these shared operations:

- [packages/treetime/src/commands/timetree/optimization/reroot.rs](../../packages/treetime/src/commands/timetree/optimization/reroot.rs) imports `ClockParams`, `estimate_clock_model_with_reroot_policy()`, `BranchPointOptimizationParams`, `RerootChanges`, `RerootParams`
- [packages/treetime/src/commands/timetree/partition_ops.rs#L14-L24](../../packages/treetime/src/commands/timetree/partition_ops.rs#L14-L24) defines `PartitionRerootOps` trait with `apply_reroot()` -- depends on `RerootChanges` from `commands/clock/`
- [packages/treetime/src/commands/timetree/run.rs](../../packages/treetime/src/commands/timetree/run.rs) imports `RerootParams`, `estimate_clock_model_with_reroot_policy()`, `ClockParams`
- [packages/treetime/src/commands/timetree/refinement.rs#L97-L105](../../packages/treetime/src/commands/timetree/refinement.rs#L97-L105) calls `estimate_clock_model_with_reroot()`

Similarly, `optimize` imports `merge_shared_mutation_branches` from `commands/prune/`.

## Independent reimplementations

[packages/treetime/src/commands/ancestral/\_\_tests\_\_/test_marginal_root_invariance_prop.rs#L75-L117](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_root_invariance_prop.rs#L75-L117) defines `reroot_at_internal_node()` -- a standalone reroot implementation that duplicates the edge-inversion logic from `apply_reroot()` in `clock/reroot.rs`.

[packages/treetime/src/commands/ancestral/\_\_tests\_\_/test_marginal_root_invariance_prop.rs#L124-L157](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_root_invariance_prop.rs#L124-L157) defines `bypass_degree2_node()` -- functionally equivalent to `remove_node_if_trivial()` in `clock/reroot.rs`, with different error handling and API shape.

These test-local reimplementations cannot use the `clock/reroot.rs` versions because:

- `reroot_in_place()` requires `ClockNode + ClockEdge` bounds (the test graph uses `NodeAncestral`/`EdgeAncestral`)
- The generic graph operations are entangled with clock-specific orchestration

## Architectural inconsistency

A shared `topology_cleanup/` module already exists at [packages/treetime/src/representation/algo/topology_cleanup/](../../packages/treetime/src/representation/algo/topology_cleanup/) with `collapse_edge()` [packages/treetime/src/representation/algo/topology_cleanup/collapse.rs#L33](../../packages/treetime/src/representation/algo/topology_cleanup/collapse.rs#L33). This module is the designated home for cross-command graph operations. The reroot operations have not been extracted to follow the same pattern.

## Layering

The current reroot call chain has two layers:

1. **Graph-topology layer** (generic): edge splitting, path inversion, degree-2 node removal, edge merging. No clock or partition awareness. Currently in `clock/reroot.rs`.
2. **Clock-aware layer**: `estimate_clock_model_with_reroot_policy()` [packages/treetime/src/commands/clock/clock_regression.rs#L161-L220](../../packages/treetime/src/commands/clock/clock_regression.rs#L161-L220) runs clock regression, calls `reroot_in_place()`, extracts the clock model.
3. **Partition-aware layer**: `reroot_tree()` [packages/treetime/src/commands/timetree/optimization/reroot.rs#L20-L96](../../packages/treetime/src/commands/timetree/optimization/reroot.rs#L20-L96) wraps the clock layer, then propagates topology changes to partitions via `PartitionRerootOps::apply_reroot()`.

The partition-aware `apply_reroot()` implementation for sparse partitions is at [packages/treetime/src/representation/partition/marginal_sparse.rs#L133](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L133), which is correctly placed in the representation layer.

The problem is layer 1: generic graph operations are embedded in layer 2's module.

## Impact

- Adding a command that needs reroot requires importing from `commands/clock/`, creating a false dependency on clock logic
- Test code that needs reroot without clock bounds must reimplement the operations
- The `bypass_degree2_node()` test reimplementation diverges from `remove_node_if_trivial()` in edge-case handling, risking test validity
- The `topology_cleanup/` extraction pattern is incomplete: `collapse_edge` was extracted, but reroot operations were not

## Proposed solutions

### S1. Extract graph-topology operations to `topology_cleanup/`

Move the generic operations (edge splitting, path inversion, degree-2 node removal) and their data types to `representation/algo/topology_cleanup/reroot.rs`, alongside `collapse.rs`. Keep `reroot_in_place()` in `clock/reroot.rs` as a thin clock-specific wrapper that calls `find_best_root()` then delegates to the shared topology functions.

Advantages:

- Consistent with the existing `collapse_edge` extraction
- `ancestral` tests can import shared operations directly
- Graph-topology types (`EdgeSplitInfo`, `EdgeMergeInfo`, `RerootChanges`) live where they belong

### S2. Extract to `treetime-graph` crate

Move the generic graph operations to the `treetime-graph` crate, where `invert_edge()` already lives. The reroot data types and operations are pure graph-topology concerns with no domain dependencies.

Advantages:

- Operations live at the same abstraction level as `invert_edge()`
- Eliminates cross-crate coupling for graph operations
- Any crate in the workspace can use them without depending on `treetime`

Disadvantage:

- `treetime-graph` is currently a thin graph data structure crate; adding algorithms broadens its scope

### S3. Shared `commands/shared/` module

Create `commands/shared/reroot.rs` for the extracted operations. This is the approach suggested in [L-core-duplicate-write-graph.md](L-core-duplicate-write-graph.md) for `write_graph()`.

Advantages:

- Low-friction extraction (stays within the `treetime` crate)
- Groups all cross-command shared code in one place

Disadvantage:

- Graph topology operations are not command-level concerns; `representation/algo/` is a better fit

## Related

- [L-core-duplicate-write-graph.md](L-core-duplicate-write-graph.md) -- same pattern: `write_graph()` duplicated across `ancestral`, `optimize`, `prune`
- [L-topology-cleanup-move-merge-shared-mutations.md](L-topology-cleanup-move-merge-shared-mutations.md) -- same pattern: `merge_shared_mutation_branches()` in `prune/`, imported by `optimize/`
- [packages/treetime/src/representation/algo/topology_cleanup/](../../packages/treetime/src/representation/algo/topology_cleanup/) -- target module for graph operation extraction (already has `collapse_edge`)
