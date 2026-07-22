# Mugration confidence rows are copied for output

`PartitionMarginalDiscrete::get_confidence()` copies `node.profile.dis.row(0)` into a new `Array1<f64>`. Mugration node-data and tree-output projection call it independently, and `MugrationConfidenceOutput::new()` eagerly copies every row into a second output structure retained beside the partition.

The posterior matrix already owns these rows for the lifetime of `MugrationGraphData`. Read-only output should borrow them.

## Required behavior

- Return `Option<ArrayView1<'_, f64>>` from `get_confidence()`.
- Accept `ArrayView1<'_, f64>` in confidence-map and entropy calculations.
- Derive both values from one borrowed row per node.
- Remove the eagerly duplicated `MugrationConfidenceOutput` profile matrix. Render confidence CSV directly from graph node names and partition row views at write time.
- Keep owned arrays only when an output value must outlive the partition.

## Validation

- Exact node-data, Auspice, and confidence-CSV equivalence.
- Owned, row, and non-contiguous view unit cases for consumers.
- Allocation regression coverage proving projection does not copy confidence rows.

## Related tickets

- [kb/tickets/mugration-borrow-confidence-rows.md](../tickets/mugration-borrow-confidence-rows.md)
