# ~~Dummy GTR initialization pattern across commands~~ (RESOLVED)

The dummy-then-replace pattern is eliminated for all commands via the `PartitionFitch` typestate lifecycle:

- **Sparse** (`ancestral`, `optimize`, `prune`, `timetree`): `PartitionFitch::compress` -> `resolve_gtr`/`infer_gtr` -> `into_marginal_sparse(gtr, &graph)` constructs the `PartitionMarginalSparse` with the real GTR from the start. No dummy needed because Fitch compression does not require a GTR model.
- **Dense `--model=infer`** (`ancestral`, `optimize`, `timetree`): `PartitionFitch::compress` -> `infer_gtr` -> `into_marginal_dense(gtr)` infers GTR from Fitch counts before constructing the dense partition. No JC69 placeholder or double marginal pass needed.
- **Dense named models** (`ancestral`, `optimize`, `timetree`): `PartitionMarginalDense::new(index, gtr, alphabet, length)` constructs with the real GTR directly. No dummy needed because named models are known upfront.

## Previous state

Three commands (`ancestral`, `optimize`, `prune`) created partition structs with a dummy JC69 GTR model, populated the partition with sequence data, then replaced the GTR with the real model. The code was annotated with FIXME comments. The root cause was a circular dependency: partitions needed a GTR at construction, but GTR inference needed populated partitions.

## Related

- [Dense partitions lack Fitch compression](H-dense-with-fitch-compression.md) - the Fitch-based GTR inference that enabled this resolution
