# Owned ndarray signatures force projection and propagation copies

> [!WARNING]
> **Needs review.** The confidence half of this issue is partly stale. Each confidence row is now copied once per node, not twice: `fn gather_reconstruction_maps()` calls `get_confidence()` once per node [packages/treetime/src/mugration/pipeline.rs#L189-L195](../../packages/treetime/src/mugration/pipeline.rs#L189-L195) and stores the owned rows in `MugrationOutput::confidences` [packages/treetime/src/mugration/pipeline.rs#L171](../../packages/treetime/src/mugration/pipeline.rs#L171). The output writers read those stored rows. A borrowed return from `get_confidence()` only removes the copy if the output keeps the partition alive instead of owning the rows, which is the scope of [N-mugration-confidence-rows-copied-for-output.md](N-mugration-confidence-rows-copied-for-output.md).

Read-only APIs accept owned ndarray references or return owned arrays even when their data already lives in a longer-lived matrix. Callers must copy confidence rows and transposed transition matrices at per-node, per-edge, or per-site frequency.

## Evidence

- `PartitionMarginalDiscrete::get_confidence()` [packages/treetime/src/partition/marginal/discrete/partition.rs#L78-L85](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L78-L85) returns `Array1<f64>` by copying `node.profile.dis.row(0)`.
- `fn build_confidence_map()` and `fn compute_entropy()` [packages/app-output/src/mugration_tree_output.rs#L167-L180](../../packages/app-output/src/mugration_tree_output.rs#L167-L180) accept `&Array1<f64>`, rejecting the row view naturally produced by the partition. Callers: [packages/app-output/src/mugration_tree_output.rs#L122-L126](../../packages/app-output/src/mugration_tree_output.rs#L122-L126) and [packages/app-output/src/augur_node_data_mugration.rs#L115-L121](../../packages/app-output/src/augur_node_data_mugration.rs#L115-L121).
- `fn propagate_raw()` [packages/treetime/src/partition/marginal/sparse/message.rs#L105-L109](../../packages/treetime/src/partition/marginal/sparse/message.rs#L105-L109) accepts `&Array2<f64>`. Its backward caller converts the zero-copy transpose view with `.t().to_owned()` [packages/treetime/src/partition/marginal/sparse/backward.rs#L149-L153](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L149-L153).
- `fn propagate_raw_per_site()` copies the transposed matrix once per call and once per variable site [packages/treetime/src/partition/marginal/sparse/message.rs#L157-L182](../../packages/treetime/src/partition/marginal/sparse/message.rs#L157-L182).

## Options

- **Concrete view types:** accept `ArrayView1`, `ArrayView2`, `ArrayViewMut1`, and `ArrayViewMut2`. Lifetimes and mutability remain explicit, and common row, slice, and transpose callers need no copy.
- **Generic storage bounds:** accept `ArrayBase<S, D>`. This supports owned arrays and views through one signature, but exposes more generic parameters where the function needs only a borrowed view.

## Recommendation

Use concrete view types at borrowing boundaries. Return `Option<ArrayView1<'_, f64>>` from `get_confidence()`, bind it once per node during projection, and derive both confidence and entropy from it. Accept `ArrayView2` in `propagate_raw()` and the other read-only propagation boundaries where storage ownership is irrelevant, and pass transpose and sliced views directly. Keep owned return values only where the callee constructs independent data.

## Required properties

- Owned, contiguous-view, non-contiguous sliced-view, row-view, and transposed-view inputs produce identical values.
- Projection obtains one borrowed confidence row per node and derives both confidence and entropy from it.
- Backward propagation does not allocate a copied transition matrix solely to transpose it.

## Validation

- Unit cases for owned, row-view, transposed-view, and non-contiguous sliced-view inputs.
- Exact output equivalence and allocation regression tests.

## Related

- [N-ancestral-marginal-array-kernels-allocate.md](N-ancestral-marginal-array-kernels-allocate.md)
- [N-mugration-confidence-rows-copied-for-output.md](N-mugration-confidence-rows-copied-for-output.md): overlaps on the borrowed `get_confidence()` return and the confidence and entropy consumers
