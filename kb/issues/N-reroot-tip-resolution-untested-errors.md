# resolve_tip_keys error and ambiguity paths untested

> [!IMPORTANT]
> **Investigation required.** Expected behavior for empty input, duplicate names, and internal-node matches must be established before complete tests or a shared resolver can be specified.

`fn resolve_tip_keys()` [packages/treetime/src/optimize/pipeline.rs#L251](../../packages/treetime/src/optimize/pipeline.rs#L251) is only tested with valid leaf names. Missing coverage includes a missing tip name, an empty tip list, duplicate names in the tree, and internal-node name matches. The lookup is duplicated by `fn find_node_key_by_name()` [packages/treetime/src/clock/reroot.rs#L270](../../packages/treetime/src/clock/reroot.rs#L270), tracked in [kb/issues/N-reroot-duplicated-tip-name-resolution.md](N-reroot-duplicated-tip-name-resolution.md).

## Related KB items

- [kb/issues/N-reroot-duplicated-tip-name-resolution.md](N-reroot-duplicated-tip-name-resolution.md)
- [kb/proposals/reroot-generic-scoring-architecture.md](../proposals/reroot-generic-scoring-architecture.md)
- [kb/proposals/optimize-reroot-support.md](../proposals/optimize-reroot-support.md)
