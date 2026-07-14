# Owned ndarray signatures force projection and propagation copies

Read-only APIs accept owned ndarray references or return owned arrays even when their data already lives in a longer-lived matrix. Callers must copy confidence rows and transposed transition matrices at per-node, per-edge, or per-site frequency.

## Evidence

- `MarginalDiscretePartition::get_confidence()` [packages/treetime/src/partition/marginal_discrete.rs#L99](../../packages/treetime/src/partition/marginal_discrete.rs#L99) returns `Array1<f64>` by copying `node.profile.dis.row(0)`.
- `fn build_ir_mugration()` [packages/treetime/src/commands/shared/ir_projection.rs#L126](../../packages/treetime/src/commands/shared/ir_projection.rs#L126) calls that accessor twice per node, once for the confidence map and once for entropy.
- `fn build_confidence_map()` and `fn compute_entropy()` [packages/treetime/src/commands/shared/ir_projection.rs#L197](../../packages/treetime/src/commands/shared/ir_projection.rs#L197) accept `&Array1<f64>`, rejecting the row view naturally produced by the partition.
- `fn propagate_raw()` [packages/treetime/src/partition/marginal_helpers.rs#L90](../../packages/treetime/src/partition/marginal_helpers.rs#L90) accepts `&Array2<f64>`. Its callers in [packages/treetime/src/partition/marginal_passes.rs#L399](../../packages/treetime/src/partition/marginal_passes.rs#L399) and [packages/treetime/src/partition/marginal_helpers.rs#L138](../../packages/treetime/src/partition/marginal_helpers.rs#L138) convert zero-copy transpose views with `.t().to_owned()`.

## Options

- **Concrete view types:** accept `ArrayView1`, `ArrayView2`, `ArrayViewMut1`, and `ArrayViewMut2`. Lifetimes and mutability remain explicit, and common row, slice, and transpose callers need no copy.
- **Generic storage bounds:** accept `ArrayBase<S, D>`. This supports owned arrays and views through one signature, but exposes more generic parameters where the function needs only a borrowed view.

## Recommendation

Use concrete view types at borrowing boundaries. Return `Option<ArrayView1<'_, f64>>` from `get_confidence()`, bind it once during projection, and pass transpose views directly into propagation. Keep owned return values only where the callee constructs independent data.

## Required properties

- Owned, contiguous-view, non-contiguous sliced-view, row-view, and transposed-view inputs produce identical values.
- Projection obtains one borrowed confidence row per node and derives both confidence and entropy from it.
- Backward propagation does not allocate a copied transition matrix solely to transpose it.

## Related issues

- [N-ancestral-marginal-array-kernels-allocate.md](N-ancestral-marginal-array-kernels-allocate.md)
