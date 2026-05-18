# Polytomy resolution test improvements

Test coverage gaps and quality issues in polytomy resolution tests.

v1 tests: [`packages/treetime/src/timetree/optimization/__tests__/test_polytomy.rs`](../../packages/treetime/src/timetree/optimization/__tests__/test_polytomy.rs)

## Items

### Point distributions make 4 topology tests non-exercising of cost function

`create_polytomy_tree()` uses `Distribution::point(0.1, 1.0)` on edges. Point distributions error for any `t != 0.1`, so all `eval()` calls fall through to `.unwrap_or(1e-10)`. The cost function output is independent of distribution shape. Tests pass because `threshold = -1000.0` accepts any negative gain. The 4 affected tests verify topology operations correctly but do not exercise the cost function. `create_polytomy_tree_with_realistic_distributions()` (used by 2 slope tests) shows the correct approach.

### Penalty formula tested directionally only

Tests verify that higher `zero_branch_slope` suppresses merges and zero slope allows them, but no test verifies the linear formula `penalty = zero_branch_slope * dt`. A change to `slope * dt^2` would pass all existing tests.

### Guard path untested

`compute_merge_gain()` returns `None` when `parent_time >= child_min_time`. No test exercises this guard. A test with equal or reversed parent/child times would cover this.

### Edge `time_length` values not verified

`merge_children()` creates edges with `time_length` but no test verifies the arithmetic (`child_time - new_node_time` for child edges, `new_node_time - parent_time` for parent edge).

### No structural integrity checks after resolution

Tests check degree counts but not leaf reachability, edge count consistency, or DAG validity after topology-modifying operations.

### Tautological time assertion

`test_resolve_polytomies_new_node_has_correct_time` asserts new node time is in `(parent_time, min_child_time)`, which is the exact domain enforced by Brent optimizer bounds. With point distributions (above), the assertion is trivially satisfied. With realistic distributions, the assertion tests a necessary but not sufficient postcondition.
