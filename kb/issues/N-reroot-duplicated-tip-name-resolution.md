# Tip-name resolution duplicated between optimize and clock reroot

Two reroot call sites resolve the tip names of `--reroot-tips` to node keys with the same code. Optimize's `fn resolve_tip_keys()` [packages/treetime/src/optimize/pipeline.rs#L278](../../packages/treetime/src/optimize/pipeline.rs#L278) and clock's `fn find_tip_group_root()` [packages/treetime/src/clock/reroot.rs#L283](../../packages/treetime/src/clock/reroot.rs#L283) both scan the name map for the first node whose name equals the tip and fail with `Reroot tip not found: <tip>` otherwise:

```rust
names
  .iter()
  .find(|(_, name)| name.as_deref() == Some(tip.as_str()))
  .map(|(key, _)| *key)
  .ok_or_else(|| input_error(format!("Reroot tip not found: {tip}")))
```

## Duplicate names

The name map is ordered by node key, and node keys follow the input order of the Newick tree, so a duplicate tip name resolves to its first node in input order. This is the rule of [kb/decisions/duplicate-names-warned-ids-from-input-order.md](../decisions/duplicate-names-warned-ids-from-input-order.md). Neither call site warns. Under that decision, a warning is required only when the duplicate is already known or cheap to find at this point, and the lookup must not add a duplicate search that every run pays for.

## Required behavior

Extract one shared tip lookup that keeps the first-match rule, and use it from both reroot paths. The migration of clock rerooting onto the generic reroot module is the natural consolidation point.

## Locations

- `fn resolve_tip_keys()` [packages/treetime/src/optimize/pipeline.rs#L278](../../packages/treetime/src/optimize/pipeline.rs#L278)
- `fn find_tip_group_root()` [packages/treetime/src/clock/reroot.rs#L283](../../packages/treetime/src/clock/reroot.rs#L283)
- [kb/proposals/reroot-generic-scoring-architecture.md](../proposals/reroot-generic-scoring-architecture.md)
- [kb/issues/N-reroot-tip-resolution-untested-errors.md](N-reroot-tip-resolution-untested-errors.md)
