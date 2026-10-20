# `infer_dense()` stub always returns false

> [!IMPORTANT]
> **Decision required.** The two fix options under "Fix" lead to different code: implement a selection heuristic, or remove the shared selector and state the default where each pipeline resolves its representation. [kb/decisions/sequence-representation-dense-sparse.md](../decisions/sequence-representation-dense-sparse.md) describes the planned criterion (expected mutations per branch relative to sequence length) but no threshold. A third option, until sparse and dense agree, is to select dense for small inputs, where its memory cost is low.

`fn infer_dense()` [packages/treetime/src/partition/algo/infer_dense.rs#L1-L7](../../packages/treetime/src/partition/algo/infer_dense.rs#L1-L7) is the shared dense-vs-sparse selector, but the function always returns `false`. `fn Representation::resolve()` uses it as the default when `--dense` is omitted, and the `ancestral`, `optimize`, and `timetree` pipelines call `Representation::resolve()`.

## Impact

Omitting `--dense` always selects sparse mode. On long-branch datasets where dense representation would be more appropriate, users must explicitly pass `--dense=true`. Sparse mode is the default of every run, so the differences between sparse and dense results in [H-timetree-sparse-dense-node-times-differ.md](H-timetree-sparse-dense-node-times-differ.md), [H-ancestral-sparse-dense-internal-residues-diverge.md](H-ancestral-sparse-dense-internal-residues-diverge.md), and [H-sparse-marginal-likelihood-depends-on-branch-of-unknown-leaf.md](H-sparse-marginal-likelihood-depends-on-branch-of-unknown-leaf.md) reach every run that omits `--dense`.

## Root cause

The heuristic for automatic selection has not been implemented. The stub was created as a placeholder during initial command wiring.

## Fix

Implement a heuristic based on branch lengths, sequence length, or both. Alternatively, remove the shared abstraction and make the default explicit at each command until a real heuristic exists.

## Locations

- Stub: `fn infer_dense()` [packages/treetime/src/partition/algo/infer_dense.rs#L1-L7](../../packages/treetime/src/partition/algo/infer_dense.rs#L1-L7)
- Selector: `fn Representation::resolve()` [packages/treetime/src/partition/create.rs#L64-L70](../../packages/treetime/src/partition/create.rs#L64-L70)
- Consumers:
  - ancestral, nucleotide [packages/treetime/src/ancestral/plan.rs#L49](../../packages/treetime/src/ancestral/plan.rs#L49)
  - ancestral, amino acid [packages/treetime/src/ancestral/aa.rs#L36](../../packages/treetime/src/ancestral/aa.rs#L36)
  - optimize [packages/treetime/src/optimize/pipeline.rs#L66](../../packages/treetime/src/optimize/pipeline.rs#L66)
  - timetree [packages/treetime/src/timetree/pipeline.rs#L328](../../packages/treetime/src/timetree/pipeline.rs#L328)

## Related decisions

- [kb/decisions/sequence-representation-dense-sparse.md](../decisions/sequence-representation-dense-sparse.md)
