# Mugration confidence rows are copied for output

`PartitionMarginalDiscrete::get_confidence()` copies `node.profile.dis.row(0)` into a new `Array1<f64>`. Mugration node-data and tree-output projection call it independently, and `MugrationConfidenceOutput::new()` eagerly copies every row into a second output structure retained beside the partition.

The posterior matrix already owns these rows for the lifetime of `MugrationGraphData`. Read-only output should borrow them.

## Required behavior

- Return `Option<ArrayView1<'_, f64>>` from `get_confidence()`.
- Accept `ArrayView1<'_, f64>` in confidence-map and entropy calculations.
- Derive both values from one borrowed row per node.
- Remove the eagerly duplicated `MugrationConfidenceOutput` profile matrix. Render confidence CSV directly from graph node names and partition row views at write time.
- Update every `get_confidence()` caller without adding ownership adapters.
- Keep owned arrays only when an output value must outlive the partition.

## Validation

- Exact node-data, Auspice, and confidence-CSV equivalence.
- Owned, row, and non-contiguous view unit cases for consumers.
- Allocation regression coverage proving projection does not copy confidence rows.

## Locations

- `fn PartitionMarginalDiscrete::get_confidence()` [packages/treetime/src/partition/marginal/discrete/partition.rs#L78-L85](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L78-L85)
- `MugrationConfidenceOutput::new()` call [packages/app-output/src/mugration_result.rs#L44](../../packages/app-output/src/mugration_result.rs#L44), `struct MugrationConfidenceOutput` [packages/app-output/src/mugration_result.rs#L55](../../packages/app-output/src/mugration_result.rs#L55)
- `fn build_confidence_map()` and `fn compute_entropy()` [packages/app-output/src/mugration_tree_output.rs#L167-L180](../../packages/app-output/src/mugration_tree_output.rs#L167-L180)

## Related issues

- [M-io-auspice-entropy-perturbs-shannon-definition.md](M-io-auspice-entropy-perturbs-shannon-definition.md): changes the formula of `compute_entropy()`
