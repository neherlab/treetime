# Make parallel clock output order deterministic

Remove scheduling-dependent clock CSV row order while preserving parallel computation.

## Required changes

- Construct regression rows in slots keyed by the command's selected topology/reference-order index.
- Serialize rows in index order, independent of worker completion order.
- Keep CSV output stable across worker counts and repeated runs.

## Validation

- Byte comparison of repeated one-, two-, eight-, and sixteen-worker CSV outputs.
- Duplicate-slot runs of the identical binary.
- Semantic comparison against the selected topology/reference ordering.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/M-clock-parallel-output-order-nondeterministic.md](../issues/M-clock-parallel-output-order-nondeterministic.md)
