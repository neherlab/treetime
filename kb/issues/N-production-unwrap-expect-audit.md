# Production unwrap/expect/assert audit

## Summary

Production `unwrap()`/`expect()`/`assert!()` calls remain that can panic instead of returning contextual errors. The former `InDel::new()` assertions have been replaced by a fallible constructor and are no longer part of this audit.

## Production unwrap/expect/assert instances (~75+)

| File                                                 | Instances | Pattern                                                       |
| :--------------------------------------------------- | :-------- | :------------------------------------------------------------ |
| `partition/marginal_passes.rs`                       | 10        | `.unwrap()` / `.expect()` on node/edge lookup                 |
| `partition/marginal_dense.rs`                        | 8         | `.unwrap()` / `.expect()` on node/edge lookup                 |
| `clock/reroot.rs`                                    | 6         | `.expect("Edge not found")`                                   |
| `timetree/optimization/polytomy.rs`                  | 5         | `.expect("Node must exist")`                                  |
| `partition/traits.rs`                                | 4         | `.expect()` on `node()`, `edge()`, `node_mut()`, `edge_mut()` |
| `packages/treetime-grid/src/grid.rs`                 | 5         | `T::from(...).unwrap()` numeric conversions                   |
| `packages/treetime-grid/src/grid_fn.rs`              | 3         | `.unwrap()` on numeric conversions                            |
| `clock/clock_regression.rs`                          | 1         | `.expect()` on branch_length (2 traversal unwraps fixed)      |
| `clock/date_constraints.rs`                          | 2         | `.unwrap()` on name/dates                                     |
| `gtr/gtr.rs`                                         | 3         | `assert!()` / `assert_eq!()` in constructors                  |
| `gtr/get_gtr.rs`                                     | 1         | `.expect()` on JSON serialization                             |
| `packages/treetime-utils/src/datetime/datetime.rs`   | 3         | `.unwrap()` / `.expect()`                                     |
| `packages/treetime-utils/src/datetime/date_range.rs` | 4         | `.unwrap()` in `from_ymd()`                                   |
| `seq/mutation.rs`                                    | 2         | `.unwrap()` on byte access                                    |
| `seq/div.rs`                                         | 2         | `.unwrap()` on graph lookups                                  |
| `ancestral/fitch.rs`                                 | 3         | `.unwrap()` on node operations                                |
| `packages/treetime-cli/src/cli/verbosity.rs`         | 1         | `.unwrap()` on parse                                          |
| `cli/rtt_chart.rs`                                   | 1         | `assert!(!results.is_empty())`                                |

## Impact

- Graph lookup unwraps panic if tree structure is inconsistent (e.g., after failed topology operations)
- Numeric conversion unwraps panic on out-of-range values
- Production users see panic backtraces instead of actionable error messages

## Fix

Prioritize remaining calls by verified reachability: graph traversal unwraps and numeric conversion unwraps require contextual error propagation when normal input can reach them.

## Related tickets

- [kb/tickets/safety-audit-production-unwrap-expect-assert-calls.md](../tickets/safety-audit-production-unwrap-expect-assert-calls.md)
