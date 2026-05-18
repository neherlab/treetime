# Audit ~75 production unwrap/expect/assert calls and convert to Result

## Description

Approximately 75 production `unwrap()`/`expect()`/`assert!()` calls panic instead of returning errors. Graph lookup unwraps and numeric conversion unwraps are highest priority since they are reachable from normal input paths.

## Instances by file

| File                                            | Instances | Pattern                                                       |
| :---------------------------------------------- | :-------- | :------------------------------------------------------------ |
| `representation/partition/marginal_passes.rs`   | 10        | `.unwrap()` / `.expect()` on node/edge lookup                 |
| `representation/partition/marginal_dense.rs`    | 8         | `.unwrap()` / `.expect()` on node/edge lookup                 |
| `commands/clock/reroot.rs`                      | 6         | `.expect("Edge not found")`                                   |
| `commands/timetree/optimization/polytomy.rs`    | 5         | `.expect("Node must exist")`                                  |
| `representation/partition/traits.rs`            | 4         | `.expect()` on `node()`, `edge()`, `node_mut()`, `edge_mut()` |
| `treetime-grid/src/grid.rs`                     | 5         | `T::from(...).unwrap()` numeric conversions                   |
| `treetime-grid/src/grid_fn.rs`                  | 3         | `.unwrap()` on numeric conversions                            |
| `commands/clock/clock_regression.rs`            | 3         | `.expect()` / `.unwrap()`                                     |
| `commands/clock/date_constraints.rs`            | 2         | `.unwrap()` on name/dates                                     |
| `gtr/gtr.rs`                                    | 3         | `assert!()` / `assert_eq!()` in constructors                  |
| `gtr/get_gtr.rs`                                | 1         | `.expect()` on JSON serialization                             |
| `treetime-utils/src/datetime/datetime.rs`       | 3         | `.unwrap()` / `.expect()`                                     |
| `treetime-utils/src/datetime/date_range.rs`     | 4         | `.unwrap()` in `from_ymd()`                                   |
| `seq/mutation.rs`                               | 2         | `.unwrap()` on byte access                                    |
| `seq/div.rs`                                    | 2         | `.unwrap()` on graph lookups                                  |
| `seq/indel.rs`                                  | 1         | `assert!()` on range validity                                 |
| `commands/ancestral/fitch.rs`                   | 3         | `.unwrap()` on node operations                                |
| `coalescent/contributions.rs` | 1         | `.unwrap()` on Result in parallel traversal                   |
| `commands/timetree/inference/backward_pass.rs`  | 1         | `.unwrap()` on Result                                         |
| `commands/timetree/inference/forward_pass.rs`   | 1         | `.unwrap()` on Result                                         |
| `treetime-cli/src/cli/verbosity.rs`             | 1         | `.unwrap()` on parse                                          |
| `cli/rtt_chart.rs`                              | 1         | `assert!(!results.is_empty())`                                |

The inference traversal unwraps are tracked separately in `M-timetree-inference-unwrap-in-traversals.md`.

## Impact

- Graph lookup unwraps panic if tree structure is inconsistent (e.g., after failed topology operations)
- Numeric conversion unwraps panic on out-of-range values
- Production users see panic backtraces instead of actionable error messages

## Fix

Prioritize by reachability: graph traversal unwraps and numeric conversion unwraps first (reachable from normal input paths). Replace with proper `Result` propagation and actionable error messages.

## Related issues

- Source: [N-production-unwrap-expect-audit.md](../issues/N-production-unwrap-expect-audit.md) -- delete after full resolution
- [N-error-suppression-unwrap-or-defaults.md](../issues/N-error-suppression-unwrap-or-defaults.md) -- related error suppression patterns
