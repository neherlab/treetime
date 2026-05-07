# Rebuild Alphabet serde skip fields after deserialization

Two related issues in the alphabet module.

## Instances

### serde skip fields not rebuilt

v1: [`packages/treetime/src/alphabet/alphabet.rs#L48-L56`](../../packages/treetime/src/alphabet/alphabet.rs#L48-L56)

Four `#[serde(skip)]` fields with no post-deserialization rebuild. If `struct Alphabet` is deserialized, lookup tables are empty/default.

### Missing empty-canonical rejection

v1: [`packages/treetime/src/alphabet/alphabet_config.rs#L72`](../../packages/treetime/src/alphabet/alphabet_config.rs#L72)

`fn validate()` does not reject an empty canonical character set. An empty alphabet produces division-by-zero in profile construction.

## Related issues

- Source: [N-code-quality-conventions.md](../issues/N-code-quality-conventions.md) -- delete after full resolution
