# Fitch passes leave the partition without nodes on error

`fn fitch_backward()` and `fn fitch_forward()` move the node map out of the partition with `mem::take(&mut partition.nodes)` before they run a fallible graph pass. When a visitor returns an error, the `?` returns before the outputs are written back, so `partition.nodes` stays empty while `partition.edges` keeps its previous value:

- `fn fitch_backward()` [packages/treetime/src/partition/fitch/passes.rs#L91-L104](../../packages/treetime/src/partition/fitch/passes.rs#L91-L104)
- `fn fitch_forward()` [packages/treetime/src/partition/fitch/passes.rs#L172-L188](../../packages/treetime/src/partition/fitch/passes.rs#L172-L188)

Both passes run from `fn compress_sequences()` [packages/treetime/src/partition/fitch/passes.rs#L43-L52](../../packages/treetime/src/partition/fitch/passes.rs#L43-L52), which borrows the partition mutably, so a caller can keep the emptied partition after the error.

The marginal, branch-optimization, and timetree branch-distribution passes do not have this defect: they compute complete result maps and publish them only after every worker succeeds.

## Impact

The command returns an error but leaves a partition whose node map does not match its edge map. Retrying or inspecting the partition after the failure is unsafe.

## Potential solutions

- O1. Run the pass on borrowed nodes (`map_backward`, `map_forward`) and replace `partition.nodes` and `partition.edges` only after success. The owned variants exist to avoid copying node sequences (commit `b7cfe761`), so this option must show that the borrowed pass does not reintroduce those copies.
- O2. Keep the owned pass and give it back its unconsumed inputs on failure, so the caller can restore the node map. Nodes already consumed by successful visitors must be recoverable as well, which the current slot design does not provide.

## Recommendation

Commit nodes and edges together only after the pass succeeds. The choice between O1 and O2 depends on the copy cost of O1, which is not measured.

## Validation

- Inject a visitor failure at the first, a middle, and the last node under one and several Rayon worker counts, and compare the complete partition with its state before the call.
- Assert that successful results are identical across worker counts and repeated runs.

## Related issues

- [N-ancestral-parallel-sparse-leaf-error-atomicity-unverified.md](N-ancestral-parallel-sparse-leaf-error-atomicity-unverified.md)
