# Borrow mugration confidence rows during output

Use ndarray row views throughout mugration output and remove the duplicated confidence matrix retained in `MugrationGraphData`.

## Required changes

- Return `Option<ArrayView1<'_, f64>>` from `PartitionMarginalDiscrete::get_confidence()`.
- Accept `ArrayView1<'_, f64>` in `build_confidence_map()` and exact Shannon entropy calculation.
- Bind one row view per node and derive confidence plus entropy from it.
- Replace `MugrationConfidenceOutput`'s owned rows with write-time CSV traversal over final graph names and partition views.
- Update every `get_confidence()` caller without adding ownership adapters.

## Validation

- Exact confidence CSV, Augur node-data, and Auspice output equivalence.
- Deterministic, uniform, and non-contiguous profile view cases.
- Allocation regression coverage for per-node output projection.
- Full lint and test suite.

## Related issues

Source: [kb/issues/N-mugration-confidence-rows-copied-for-output.md](../issues/N-mugration-confidence-rows-copied-for-output.md)
