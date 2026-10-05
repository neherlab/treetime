# Mugration confidence rows are copied for output

`PartitionMarginalDiscrete::get_confidence()` copies `node.profile.dis.row(0)` into a new `Array1<f64>`. The mugration pipeline calls it for every node and keeps the copies in `MugrationOutput::confidences`, beside the posterior matrix that already owns these rows for the lifetime of `MugrationGraphData`. Read-only output should borrow them.

The output layer reads the copies once: Auspice JSON, augur node data and the confidence CSV take each node's profile from `TreeTraits::profiles` of the struct of facts, and compute the confidence map and the entropy from it without a further copy.

## Required behavior

- Return `Option<ArrayView1<'_, f64>>` from `get_confidence()`.
- Accept `ArrayView1<'_, f64>` in `build_confidence_map()` and `compute_entropy()` and in `TreeTraits::profiles`.
- Update every `get_confidence()` caller without adding ownership adapters.
- Keep owned arrays only when an output value must outlive the partition.

## Validation

- Exact node-data, Auspice, and confidence-CSV equivalence.
- Owned, row, and non-contiguous view unit cases for consumers.
- Allocation regression coverage proving projection does not copy confidence rows.

## Locations

- `fn PartitionMarginalDiscrete::get_confidence()` [packages/treetime/src/partition/marginal/discrete/partition.rs#L78-L85](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L78-L85)
- `MugrationOutput::confidences` filled in [packages/treetime/src/mugration/pipeline.rs#L193](../../packages/treetime/src/mugration/pipeline.rs#L193)
- `struct TreeTraits` in [packages/app-output/src/annotated_graph.rs](../../packages/app-output/src/annotated_graph.rs)
- `fn build_confidence_map()` and `fn compute_entropy()` in [packages/app-output/src/trait_profile.rs](../../packages/app-output/src/trait_profile.rs)

## Related issues

- [M-io-auspice-entropy-perturbs-shannon-definition.md](M-io-auspice-entropy-perturbs-shannon-definition.md): changes the formula of `compute_entropy()`
