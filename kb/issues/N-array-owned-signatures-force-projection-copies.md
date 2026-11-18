# Owned ndarray signatures force transition-matrix copies

The sparse propagation function accepts owned ndarray references even when the data already lives in a longer-lived matrix. Callers must copy the transposed transition matrix once per edge.

## Evidence

- `fn propagate_raw()` [packages/treetime/src/partition/marginal/sparse/message.rs#L202](../../packages/treetime/src/partition/marginal/sparse/message.rs#L202) accepts `&Array2<f64>`. Its backward caller converts the zero-copy transpose view with `.t().to_owned()` [packages/treetime/src/partition/marginal/sparse/backward.rs#L93](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L93).

The mugration confidence rows are not part of this issue. `MugrationOutput` owns one `StatePosterior` per node with a copy of its row, and the pipeline drops the node states after the readout [packages/treetime/src/mugration/pipeline.rs#L135-L143](../../packages/treetime/src/mugration/pipeline.rs#L135-L143). A borrowed profile would have to keep every node state alive, which costs more memory than one row per node.

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


## Proposed implementation plan

1. Change `fn propagate_raw()` in `packages/treetime/src/partition/marginal/sparse/message.rs` (line 202) to take `exp_qt: ArrayView2<f64>`, and keep its body unchanged apart from the type, because it reads `exp_qt` only.
2. Update both callers: `packages/treetime/src/partition/marginal/sparse/backward.rs` line 93 passes `gtr.expQt(branch_length).t()` without `.to_owned()`, and `packages/treetime/src/partition/marginal/sparse/forward.rs` line 216 passes `gtr.expQt(branch_length).view()`. Bind the owned `expQt` result to a local first where the borrow of a temporary does not live long enough.
3. Search the other read-only propagation boundaries (`rg '&Array2<f64>' packages/treetime/src/partition/marginal`) and convert the ones where a caller builds a copy only to pass it; leave signatures whose callers already hold owned matrices.
4. Tests: an `rstest` case table in `packages/treetime/src/partition/marginal/sparse/__tests__/` that calls `propagate_raw()` with a small non-symmetric matrix as an owned view, as `.t()` of its transpose, and as a non-contiguous slice of a larger matrix, and asserts exact equality of the resulting `SparseSeqDistribution`. Oracle: the owned contiguous input; a transpose of a transpose and a slice hold the same values, so the results must be bit-identical.
5. Validate: `just lint-rs`, `just test-rs`; the sparse marginal tests must stay unchanged.
