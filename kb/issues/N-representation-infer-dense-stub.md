# `infer_dense()` stub always returns false

> [!IMPORTANT]
> **Decision required.** The two fix options under "Fix" lead to different code: implement a selection heuristic, or remove the shared selector and state the default where each pipeline resolves its representation. [kb/decisions/sequence-representation-dense-sparse.md](../decisions/sequence-representation-dense-sparse.md) describes the planned criterion (expected mutations per branch relative to sequence length) but no threshold.

`fn infer_dense()` [packages/treetime/src/partition/algo/infer_dense.rs#L1-L3](../../packages/treetime/src/partition/algo/infer_dense.rs#L1-L3) is the shared dense-vs-sparse selector, but the function always returns `false`. `fn Representation::resolve()` uses it as the default when `--dense` is omitted, and the `ancestral`, `optimize`, and `timetree` pipelines call `Representation::resolve()`.

## Impact

Omitting `--dense` always selects sparse mode. On long-branch datasets where dense representation would be more appropriate, users must explicitly pass `--dense=true`. No correctness issue: sparse mode produces valid results on all inputs.

## Root cause

The heuristic for automatic selection has not been implemented. The stub was created as a placeholder during initial command wiring.

## Fix

Implement a heuristic based on branch lengths, sequence length, or both. Alternatively, remove the shared abstraction and make the default explicit at each command until a real heuristic exists.

## Locations

- Stub: `fn infer_dense()` [packages/treetime/src/partition/algo/infer_dense.rs#L1-L3](../../packages/treetime/src/partition/algo/infer_dense.rs#L1-L3)
- Selector: `fn Representation::resolve()` [packages/treetime/src/partition/create.rs#L25-L31](../../packages/treetime/src/partition/create.rs#L25-L31)
- Consumers:
  - ancestral, nucleotide [packages/treetime/src/ancestral/plan.rs#L67](../../packages/treetime/src/ancestral/plan.rs#L67)
  - ancestral, amino acid [packages/treetime/src/ancestral/aa.rs#L37](../../packages/treetime/src/ancestral/aa.rs#L37)
  - optimize [packages/treetime/src/optimize/pipeline.rs#L61](../../packages/treetime/src/optimize/pipeline.rs#L61)
  - timetree [packages/treetime/src/timetree/pipeline.rs#L349](../../packages/treetime/src/timetree/pipeline.rs#L349)

## Related decisions

- [kb/decisions/sequence-representation-dense-sparse.md](../decisions/sequence-representation-dense-sparse.md)
