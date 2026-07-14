# Accept ndarray views and remove projection copies

Generalize read-only array APIs and stop copying confidence rows and transposed matrices.

## Required changes

- Return `Option<ArrayView1<'_, f64>>` from `MarginalDiscretePartition::get_confidence()` [packages/treetime/src/partition/marginal_discrete.rs#L99](../../packages/treetime/src/partition/marginal_discrete.rs#L99).
- In [packages/treetime/src/commands/shared/ir_projection.rs#L145](../../packages/treetime/src/commands/shared/ir_projection.rs#L145), accept `ArrayView1` in `build_confidence_map()` and `compute_entropy()`, bind one confidence view per node, and derive both outputs from it.
- Accept `ArrayView2` in `propagate_raw()` [packages/treetime/src/partition/marginal_helpers.rs#L90](../../packages/treetime/src/partition/marginal_helpers.rs#L90) and the other read-only propagation boundaries where storage ownership is irrelevant.
- Pass transpose and sliced views directly into propagation; remove `.t().to_owned()` copies used only to satisfy a concrete signature.
- Preserve owned return values only where the callee constructs new data.

## Validation

- Owned, row-view, transposed-view, and non-contiguous sliced-view cases.
- Exact output equivalence and allocation regression tests.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/N-array-owned-signatures-force-projection-copies.md](../issues/N-array-owned-signatures-force-projection-copies.md)
