# Test the common-ancestor contract directly

Add direct tests for `fn common_ancestor()` [`packages/treetime-graph/src/common_ancestor.rs#L7-L41`](../../packages/treetime-graph/src/common_ancestor.rs#L7-L41) so its graph and error contracts are not covered only through rerooting callers.

## Acceptance criteria

- Empty input returns the documented error.
- One valid node returns that node.
- Siblings return their parent.
- Nodes at different depths return their lowest common ancestor.
- Invalid node keys return actionable errors without panicking.
- Tests assert exact node identities and error classes, not only success or failure.

## Related issues

- Source: [kb/issues/N-common-ancestor-missing-direct-tests.md](../issues/N-common-ancestor-missing-direct-tests.md)
