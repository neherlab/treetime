# Implement infer_dense() heuristic or remove abstraction

`partition/algo/infer_dense.rs` exports `infer_dense()` as the shared dense-vs-sparse selector, but the function always returns `false`. Three user-facing commands use it as the default when `--dense` is omitted: `ancestral`, `optimize`, and `timetree`.

## Impact

Omitting `--dense` always selects sparse mode. On long-branch datasets where dense representation would be more appropriate, users must explicitly pass `--dense=true`. No correctness issue: sparse mode produces valid results on all inputs.

## Root cause

The heuristic for automatic selection has not been implemented. The stub was created as a placeholder during initial command wiring.

## Fix

Implement a heuristic based on branch lengths, sequence length, or both. Alternatively, remove the shared abstraction and make the default explicit at each command until a real heuristic exists.

## Locations

- Stub: `packages/treetime/src/partition/algo/infer_dense.rs`
- Consumers: `packages/treetime/src/commands/ancestral/run.rs`, `packages/treetime/src/commands/optimize/run.rs`, `packages/treetime/src/commands/timetree/initialization.rs`

## Related issues

- Source: [N-representation-infer-dense-stub.md](../issues/N-representation-infer-dense-stub.md) -- delete after full resolution
