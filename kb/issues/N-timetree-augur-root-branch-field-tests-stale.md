# Augur node-data tests still assert the root omits branch fields

Three tests in `test_augur_node_data.rs` encode the pre-`40796fd8` behaviour and fail against
the current builder. The tests are wrong: augur requires the branch fields on every node,
including the root, which is what `40796fd8` changed the builder to emit.

## Severity

Negligible — a test and documentation defect. The runtime output is the correct one; only the
assertions and a doc comment describe the superseded behaviour.

## Failing tests

[`packages/treetime/src/commands/timetree/output/__tests__/test_augur_node_data.rs`](../../packages/treetime/src/commands/timetree/output/__tests__/test_augur_node_data.rs)

| Test | Line | Stale assertion |
| --- | --- | --- |
| `test_augur_node_data_timetree_root_has_no_branch_fields` | 96-97 | `root.clock_length.is_none()` and `root.mutation_length.is_none()` |
| `test_augur_node_data_timetree_branch_length_equals_clock_length` | 71 | `assert_eq!(checked, 2)` with the comment "root has no clock_length"; the root now carries one, so three nodes are checked |
| `test_augur_node_data_timetree_mutations_mode_root_has_no_mutation_length` | 149 | `root.mutation_length.is_none()` |

## Cause

`fn build_augur_node_data_json()` previously returned `(0.0, None, None)` for a node with no
parent edge. [`40796fd8`](../../packages/treetime/src/commands/timetree/output/augur_node_data.rs#L92)
changed that to `(0.0, Some(0.0), Some(0.0))`:

```rust
None => (0.0, Some(0.0), Some(0.0)), // root node: no parent edge, so branch fields are zero
```

augur's export validates that `mutation_length` is present for every node in the node-data
JSON, the root included, so omitting the field makes the output unusable regardless of whether
the value is ever read. The tests were not updated with the builder.

## Also stale

The doc comment immediately above the changed match still describes the old behaviour
([`augur_node_data.rs#L74-L77`](../../packages/treetime/src/commands/timetree/output/augur_node_data.rs#L74-L77)):

```rust
// Per-branch fields live on the parent edge. The root has no incoming branch:
// augur emits placeholder values there, but `export v2` sets the root divergence
// to 0 regardless and never consumes the root's branch fields, so we omit them.
```

"so we omit them" is no longer true, and the reasoning it gives — that `export v2` never reads
the root's branch fields — is the reasoning that turned out to be insufficient, because
validation happens before consumption.

## Fix

Update the three assertions to expect `Some(0.0)` for the root's `clock_length` and
`mutation_length`, rename `test_augur_node_data_timetree_root_has_no_branch_fields` to reflect
that the root now carries zeroed branch fields, and rewrite the doc comment to state that the
fields are emitted as zero because augur validates their presence on every node.

## Related

- [kb/issues/N-timetree-node-data-root-branch-fields-omitted.md](N-timetree-node-data-root-branch-fields-omitted.md) — the underlying gap, resolved by the same commit
- [kb/reports/augur-node-data-json.md](../reports/augur-node-data-json.md)
