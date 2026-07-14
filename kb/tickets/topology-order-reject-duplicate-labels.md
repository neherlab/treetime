# Reject duplicate topology-order labels

Validate target-order inputs before converting them into position lookups.

## Acceptance criteria

- Reject duplicate labels in input-order, list, and reference-topology sources.
- Reject duplicate graph leaf labels.
- Reject target labels that do not identify a graph leaf.
- Preserve actionable diagnostics for missing graph leaves.
- Add parameterized negative tests for every validation class and assert that each diagnostic identifies the offending label.
- Cover duplicate positions at the first, middle, and last positions and both mean and median aggregation; test bijection validation independently of aggregation.

## Related issues

Source: [kb/issues/M-topology-order-duplicate-label-collapse.md](../issues/M-topology-order-duplicate-label-collapse.md)
