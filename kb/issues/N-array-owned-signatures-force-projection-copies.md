# Owned ndarray signatures force transition-matrix copies

The sparse propagation functions accept owned ndarray references even when the data already lives in a longer-lived matrix. Callers must copy the transposed transition matrix at per-edge and per-site frequency.

## Evidence

- `fn propagate_raw()` [packages/treetime/src/partition/marginal/sparse/message.rs#L106-L110](../../packages/treetime/src/partition/marginal/sparse/message.rs#L106-L110) accepts `&Array2<f64>`. Its backward caller converts the zero-copy transpose view with `.t().to_owned()` [packages/treetime/src/partition/marginal/sparse/backward.rs#L94-L98](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L94-L98).
- `fn propagate_raw_per_site()` copies the transposed matrix once per call and once per variable site [packages/treetime/src/partition/marginal/sparse/message.rs#L145-L169](../../packages/treetime/src/partition/marginal/sparse/message.rs#L145-L169).

The mugration confidence rows are not part of this issue. `MugrationOutput` owns one `StatePosterior` per node with a copy of its row, and the pipeline drops the node states after the readout [packages/treetime/src/mugration/pipeline.rs#L139-L148](../../packages/treetime/src/mugration/pipeline.rs#L139-L148). A borrowed profile would have to keep every node state alive, which costs more memory than one row per node.

## Options

- **Concrete view types:** accept `ArrayView2` and `ArrayViewMut2`. Lifetimes and mutability remain explicit, and transpose and slice callers need no copy.
- **Generic storage bounds:** accept `ArrayBase<S, D>`. This supports owned arrays and views through one signature, but exposes more generic parameters where the function needs only a borrowed view.

## Recommendation

Use concrete view types at borrowing boundaries. Accept `ArrayView2` in `propagate_raw()` and the other read-only propagation boundaries where storage ownership is irrelevant, and pass transpose and sliced views directly. Keep owned return values only where the callee constructs independent data.

## Required properties

- Owned, contiguous-view, non-contiguous sliced-view, and transposed-view inputs produce identical values.
- Backward propagation does not allocate a copied transition matrix solely to transpose it.

## Validation

- Unit cases for owned, transposed-view, and non-contiguous sliced-view inputs.
- Exact output equivalence and allocation regression tests.

## Related

- [N-ancestral-sparse-transition-count-allocates-per-site.md](N-ancestral-sparse-transition-count-allocates-per-site.md)
