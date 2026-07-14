# Unify tree output ordering

Make every tree-backed writer consume the same ordered topology.

## Required changes

- Apply the selected topology order once during output preparation.
- Pass the resulting graph to direct writers and every TreeIR projection.
- Remove writer-local or projection-bypassing order paths.
- Preserve the currently selected method/format availability contract; do not add or remove parsimony formats in this ticket.

## Validation

- Matrix tests for `ancestral`, `clock`, `mugration`, `optimize`, `prune`, and `timetree` across input/reference/list topology orders and every currently supported tree-backed format.
- Whole-output tests proving direct and TreeIR-backed formats contain the same ordered node topology.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/M-io-tree-backed-output-order-inconsistent.md](../issues/M-io-tree-backed-output-order-inconsistent.md)
- Related: [kb/issues/N-ancestral-auspice-json-not-produced.md](../issues/N-ancestral-auspice-json-not-produced.md)
